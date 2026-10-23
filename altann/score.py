"""Score supplied CDS models against the global-best path of one sequence.

Each model on the job's sequence and orientation is converted to oriented
zero-based half-open CDS intervals, replayed by the native core over an
N-bounded window, and reported with its delta to the global-best path on the
same window. Coordinates in the output are mapped back to the original genome.
"""
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import time
from urllib.parse import unquote

from . import __version__
from .cli import coordinates, find_core, job_from_args, run_process, thread_count

CORE_COLUMNS = ('status', 'reason', 'failed_position_1based', 'model_score', 'reference_score',
                'delta', 'window_start_1based', 'window_end_1based', 'replaced_start_1based',
                'replaced_end_1based', 'reference_gene_indices', 'identical_to_reference',
                'splice_site_jaccard')
COLUMNS = ('transcript_id', 'gene_id', 'seqid', 'strand', 'cds_start_1based', 'cds_end_1based',
           'cds_segments', 'cds_length', 'status', 'reason', 'failed_position_1based',
           'model_score', 'reference_score', 'delta', 'window_start_1based', 'window_end_1based',
           'replaced_start_1based', 'replaced_end_1based', 'reference_loci',
           'identical_to_reference', 'splice_site_jaccard')


def read_models(path):
    """Return transcript -> (gene, seqid, strand, [(start, end)]) from GFF3 or GTF CDS rows."""
    genes, cds = {}, defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9:
                raise ValueError(f'Expected nine columns in {path}: {line.strip()}')
            if 'transcript_id "' in row[8]:
                info = {k: v.strip('"') for k, _, v in
                        (item.strip().partition(' ') for item in row[8].split(';') if item.strip())}
                parents = [info['transcript_id']]
                if row[2] == 'CDS':
                    genes.setdefault(parents[0], info.get('gene_id', '.'))
            else:
                info = dict(item.split('=', 1) for item in row[8].split(';') if '=' in item)
                if row[2] in ('mRNA', 'transcript'):
                    genes[unquote(info['ID'])] = unquote(info.get('Parent', '.'))
                parents = [unquote(p) for p in info.get('Parent', '').split(',') if p]
            if row[2] != 'CDS':
                continue
            if not parents:
                raise ValueError(f'CDS without a parent transcript in {path}')
            for parent in parents:
                cds[parent].append((unquote(row[0]), row[6], int(row[3]), int(row[4])))
    models = {}
    for tid, rows in cds.items():
        if len({(seqid, strand) for seqid, strand, _, _ in rows}) != 1:
            models[tid] = (genes.get(tid, '.'), None, None, [])
            continue
        seqid, strand = rows[0][:2]
        models[tid] = (genes.get(tid, '.'), seqid, strand, sorted((s, e) for _, _, s, e in rows))
    return models


def score(args):
    started = time.monotonic()
    if args.flank < 0:
        raise ValueError('--flank must be nonnegative')
    threads = thread_count(args)
    job = job_from_args(args, require_best=False)
    if job['psauron'] or job['scores']:
        raise ValueError('score requires processed UniAnn scores (out.ps.txt and out.*.txt)')
    core = find_core(args)
    coords, segment = coordinates(job)
    chromosome, offset, length = coords['chromosome'], coords['offset'], coords['length']
    minus = job['strand'] == 'minus'
    output = args.output.resolve()
    meta_path = Path(str(output) + '.json')
    inputs = {Path(v).resolve() for k, v in job.items() if k in ('fasta', 'emissions', 'gt', 'ag', 'atg', 'stop', 'log') and v}
    inputs.add(args.models.resolve())
    if any(t in inputs for t in (output, meta_path)):
        raise ValueError('Output paths must not overwrite input files')

    counts = Counter()
    selected = []
    for tid, (gene, seqid, strand, cds) in read_models(args.models).items():
        if seqid is None:
            counts['mixed_sequence_or_strand'] += 1
        elif seqid != chromosome:
            counts['other_sequence'] += 1
        elif strand != ('-' if minus else '+'):
            counts['other_strand'] += 1
        elif cds[0][0] <= offset or cds[-1][1] > offset + length:
            counts['outside_segment'] += 1
        else:
            if minus:
                local = sorted((length - (e - offset), length - (s - offset) + 1) for s, e in cds)
            else:
                local = [(s - offset - 1, e - offset) for s, e in cds]
            selected.append((tid, gene, cds, local))
    counts['scored_models'] = len(selected)

    def genomic(position):
        if position == '.':
            return '.'
        position = int(position)
        return str(length - position + 1 + offset if minus else position + offset)

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.altann-score-', dir=output.parent) as temp:
        temp = Path(temp)
        table = temp / 'models.tsv'
        with table.open('w') as handle:
            for index, (_, _, _, local) in enumerate(selected):
                handle.write(f'm{index}\t' + ','.join(f'{s}-{e}' for s, e in local) + '\n')
        raw = temp / 'scores.tsv'
        command = [core, job['fasta'], *(job[k] for k in ('emissions', 'gt', 'ag', 'atg', 'stop')),
                   '--flank', args.flank, '--threads', threads,
                   '--score-models', table, '--output', raw]
        if job['log']:
            command.extend(['--viterbi-log', job['log']])
        run_process(command, temp / 'core.log', temp)
        with raw.open() as handle:
            rows = {row['model_id']: row for row in csv.DictReader(handle, delimiter='\t')}
        staged = temp / 'result.tsv'
        strand_name = 'minus' if minus else 'plus'
        with staged.open('w') as handle:
            handle.write('\t'.join(COLUMNS) + '\n')
            for index, (tid, gene, cds, _) in enumerate(selected):
                row = rows[f'm{index}']
                counts['status_' + row['status']] += 1
                window = sorted((genomic(row['window_start_1based']), genomic(row['window_end_1based'])),
                                key=lambda x: (x == '.', int(x) if x != '.' else 0))
                replaced = sorted((genomic(row['replaced_start_1based']), genomic(row['replaced_end_1based'])),
                                  key=lambda x: (x == '.', int(x) if x != '.' else 0))
                loci = '.' if row['reference_gene_indices'] == '.' else ','.join(
                    f'{segment}.{strand_name}.kbest.locus{i}' for i in row['reference_gene_indices'].split(','))
                values = [tid, gene, chromosome, '-' if minus else '+', str(cds[0][0]), str(cds[-1][1]),
                          str(len(cds)), str(sum(e - s + 1 for s, e in cds)), row['status'], row['reason'],
                          genomic(row['failed_position_1based']), row['model_score'], row['reference_score'],
                          row['delta'], *window, *replaced, loci, row['identical_to_reference'],
                          row['splice_site_jaccard']]
                handle.write('\t'.join(values) + '\n')
        metadata = dict(version=__version__, command=sys.argv, threads=threads, flank=args.flank,
                        reverse=args.reverse, models=str(args.models.resolve()),
                        reference_mode='traceback' if job['log'] else 'rerun',
                        inputs={k: str(v) for k, v in job.items()
                                if k in ('fasta', 'emissions', 'gt', 'ag', 'atg', 'stop', 'log') and v},
                        coordinates=coords, counts=dict(counts), seconds=time.monotonic() - started)
        staged_meta = temp / 'metadata.json'
        staged_meta.write_text(json.dumps(metadata, indent=2) + '\n')
        os.replace(staged_meta, meta_path)
        os.replace(staged, output)
    print(f'Wrote {output}', file=sys.stderr)
