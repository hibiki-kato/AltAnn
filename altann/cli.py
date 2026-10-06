"""Read existing UniAnn inputs, invoke the native decoder, and export GFF3.

Each run accepts one sequence and its matching probabilities or processed scores.
Both-strand inputs are prepared as oriented jobs before native decoding.
The native decoder reads a UniAnn traceback or reruns Viterbi when requested.
This module maps output coordinates and records provenance.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time

from . import __version__
from .convert import SEGMENT, convert

ROOT = Path(__file__).resolve().parent.parent
FILE_KEYS = ('fasta', 'psauron', 'scores', 'emissions', 'gt', 'ag', 'atg', 'stop', 'log', 'best_gff')
FILE_HELP = {
    'fasta': 'FASTA containing exactly one sequence.',
    'psauron': 'PSAURON CSV with all frame probabilities (PSAURON -a).',
    'scores': 'Site probability table: chrom pos strand type motif prob [rescaled_prob].',
    'emissions': 'Processed UniAnn emission scores (out.ps.txt).',
    'gt': 'Processed donor transition scores (out.gt.txt).',
    'ag': 'Processed acceptor transition scores (out.ag.txt).',
    'atg': 'Processed start transition scores (out.atg.txt).',
    'stop': 'Processed stop transition scores (out.stop.txt).',
    'log': 'UniAnn DP/BT log for this exact sequence and scoring model.',
    'best_gff': 'UniAnn best-path GFF matching the supplied sequence and scores.',
}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=f'AltAnn {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    d = sub.add_parser('decode', help='Decode one sequence using PSAURON/site probabilities or processed UniAnn scores.')
    d.add_argument('--input', type=Path, help='Directory containing exactly one sequence and its scores.')
    for key in FILE_KEYS:
        names = ('--best-gff', '--gff') if key == 'best_gff' else ('--' + key.replace('_', '-'),)
        if key in ('fasta', 'psauron', 'scores'):
            names += ({'fasta': '-f', 'psauron': '-p', 'scores': '-s'}[key],)
        d.add_argument(*names, type=Path, help=FILE_HELP[key])
    d.add_argument('-m', '--mult', type=float, default=2.71828182845905,
                   help='Site probability multiplier before logarithmic scaling, between 1 and 100 (default: exp(1)).')
    d.add_argument('--reverse', action='store_true', help='Input is already reverse complemented; map output coordinates back and use the negative strand.')
    d.add_argument('-a', '--all-prob', dest='both_strands', action='store_true', help='Decode both strands from the original FASTA, six-frame PSAURON CSV, and site probabilities; recompute both reference paths.')
    d.add_argument('--rerun-viterbi', action='store_true', help='Recompute the reference path from scores instead of requiring a UniAnn GFF and DP/BT log.')
    d.add_argument('--seqid', help='Chromosome identifier used in exported GFF3.')
    d.add_argument('--offset', type=int, help='Zero-based segment start in the original chromosome.')
    d.add_argument('--sequence-length', type=int, help='Full original chromosome length.')
    d.add_argument('--k', type=int, default=10, help='Retained paths per state (default: 10).')
    d.add_argument('--flank', type=int, default=1000)
    d.add_argument('--threads', default='auto', help='Native worker count or auto (available CPUs).')
    d.add_argument('--core', type=Path, help='Path to altann-core.')
    d.add_argument('--output', required=True, type=Path)
    d.add_argument('--short-cds-filter', action='store_true', help='Exclude models with <=2 exons and <=200 coding bases.')
    d.add_argument('--segment-overlap', type=int, default=4_000_000)
    d.add_argument('--segment-margin', type=int, default=20_000)
    return p


def unique(directory, patterns, label, required=False):
    found = sorted({p.resolve() for pattern in patterns for p in directory.glob(pattern) if p.is_file()})
    if len(found) > 1:
        raise ValueError(f'Ambiguous {label} in {directory}; specify its path explicitly: {found}')
    if required and not found:
        raise ValueError(f'Missing {label} in {directory}')
    return found[0] if found else None


def job_from_args(args):
    """Resolve one input set, keeping both-strand preparation separate."""
    if args.both_strands and args.reverse:
        raise ValueError('--reverse cannot be combined with --all-prob')
    if args.both_strands and args.log:
        raise ValueError('--log cannot be combined with --all-prob; both reference paths are recomputed')
    if args.rerun_viterbi and args.log:
        raise ValueError('--log cannot be combined with --rerun-viterbi')
    job = {key: getattr(args, key).resolve() if getattr(args, key) else None
           for key in FILE_KEYS}
    directory = args.input.resolve() if args.input else (job['fasta'].parent if job['fasta'] else None)
    if directory is None:
        raise ValueError('Supply --input or --fasta with its matching score files')
    if not directory.is_dir():
        raise ValueError(f'Input directory does not exist: {directory}')
    if not job['fasta']:
        job['fasta'] = unique(directory, ('*.fa', '*.fasta', '*.fna'), 'FASTA', True)
    processed = ('emissions', 'gt', 'ag', 'atg', 'stop')
    if job['psauron'] or job['scores']:
        if any(job[key] for key in processed):
            raise ValueError('PSAURON/site probabilities cannot be combined with processed score arguments')
        raw = True
    else:
        raw = not any(job[key] for key in processed) and (
            not (directory / 'out.ps.txt').is_file() or
            (args.both_strands and (directory / 'psauron_score.csv').is_file()))
    if raw:
        if args.log:
            raise ValueError('--log requires processed UniAnn score files')
        if not math.isfinite(args.mult) or not 1 <= args.mult <= 100:
            raise ValueError('--mult must be a finite number between 1 and 100')
        if not job['psauron']:
            job['psauron'] = unique(directory, ('psauron_score.csv',), 'PSAURON scores', True)
        if not job['scores']:
            job['scores'] = unique(directory, ('sites*.tsv', '*_sites.tsv', 'scores.txt'), 'site probabilities', True)
    else:
        for key, filename in [('emissions', 'out.ps.txt'), ('gt', 'out.gt.txt'),
                              ('ag', 'out.ag.txt'), ('atg', 'out.atg.txt'), ('stop', 'out.stop.txt')]:
            if not job[key]:
                job[key] = unique(directory, (filename,), key, True)
    if not job['best_gff']:
        job['best_gff'] = unique(directory, ('out.gff', 'uniann.gff', '*.uniann.gff'), 'UniAnn best GFF')
    if not args.rerun_viterbi and not args.both_strands and not raw:
        if not job['log']:
            job['log'] = unique(directory, ('out.err', 'uniann.log', '*.uniann.log', 'run.log'), 'UniAnn DP/BT log')
        if not job['best_gff'] or not job['log']:
            raise ValueError('Supply --gff and --log with a matching UniAnn DP/BT trace, '
                             'or use --rerun-viterbi to recompute the reference path')
    for key in FILE_KEYS:
        if job[key] and not job[key].is_file():
            raise ValueError(f'Missing {key} file: {job[key]}')
    job.update(strand='minus' if args.reverse else 'plus', seqid=args.seqid,
               offset=args.offset, sequence_length=args.sequence_length, mult=args.mult)
    return job


def coordinates(job):
    """Read one oriented sequence and determine its original genomic interval."""
    header, length = None, 0
    with job['fasta'].open() as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    raise ValueError('Each job requires exactly one FASTA record')
                parts = line[1:].split()
                if not parts:
                    raise ValueError('FASTA header has no sequence identifier')
                header = parts[0]
            elif line.strip():
                if header is None:
                    raise ValueError('Sequence found before FASTA header')
                length += len(''.join(line.split()))
    if not header or not length:
        raise ValueError('Empty FASTA sequence')
    match = SEGMENT.fullmatch(header)
    if match:
        seqid, offset, end, total = match.groups()
        offset, end, total = int(offset), int(end), int(total)
        if end - offset != length:
            raise ValueError(f'Segment header disagrees with FASTA length: {header}')
    else:
        seqid, offset, total = header, 0, length
    seqid = job.get('seqid') or seqid
    offset = int(job['offset']) if job.get('offset') not in (None, '') else offset
    total = int(job['sequence_length']) if job.get('sequence_length') not in (None, '') else total
    if not 0 <= offset < offset + length <= total:
        raise ValueError('Invalid offset, sequence length, or segment bounds')
    return dict(chromosome=seqid, offset=offset, length=length, total=total,
                segmented=bool(match) or offset != 0 or offset + length != total), job.get('segment') or job.get('id') or header


def run_process(command, log, cwd):
    with log.open('w') as handle:
        result = subprocess.run([str(x) for x in command], cwd=cwd, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        # Keep errors useful without printing an enormous native trace.
        with log.open('rb') as handle:
            handle.seek(max(0, log.stat().st_size - 6000))
            tail = handle.read().decode(errors='replace')
        raise ValueError(f'Command failed ({result.returncode}): {command[0]}\n{tail}')


def decode(args):
    started = time.monotonic()
    if not 1 <= args.k <= 32767 or args.flank < 0 or args.segment_overlap < 0 or args.segment_margin < 0:
        raise ValueError('K must be between 1 and 32767; flank and segment settings must be nonnegative')
    threads = (len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1) if args.threads == 'auto' else int(args.threads)
    if threads < 1:
        raise ValueError('Thread count must be positive')
    job = job_from_args(args)
    core = args.core or (ROOT / 'build' / 'altann-core')
    if not Path(core).is_file():
        installed_core = ROOT.parent.parent / 'bin' / 'altann-core'
        if args.core is None and ROOT.parent.name == 'share' and installed_core.is_file():
            core = installed_core
        else:
            core = shutil.which('altann-core') if args.core is None else None
    if not core:
        raise ValueError('altann-core was not found; build it first or specify --core')
    core = Path(core).resolve()
    output = args.output.resolve()
    targets = [output, Path(str(output)+'.tsv'), Path(str(output)+'.json')]
    inputs = {Path(v).resolve() for k, v in job.items() if k in FILE_KEYS and v}
    if any(t in inputs or (t.exists() and any(t.samefile(p) for p in inputs)) for t in targets):
        raise ValueError('Output paths must not overwrite input files')
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(version=__version__, command=sys.argv, threads=threads,
                    reverse=args.reverse, both_strands=args.both_strands,
                    reference_mode='rerun' if args.rerun_viterbi or args.both_strands or job['psauron'] else 'traceback', jobs=[])
    with tempfile.TemporaryDirectory(prefix='.altann-', dir=output.parent) as temp:
        temp = Path(temp)
        if job['psauron']:
            from .preprocess import prepare_raw
            coordinates(job)
            jobs = prepare_raw(job, temp / 'inputs', both=args.both_strands)
        elif args.both_strands:
            from .strands import prepare_both
            coordinates(job)
            jobs = prepare_both(job, temp / 'inputs')
        else:
            jobs = [job]
        database = sqlite3.connect(temp / 'merge.sqlite')
        database.execute('CREATE TABLE models (seqid TEXT, start INTEGER, finish INTEGER, strand TEXT, id TEXT, body TEXT)')
        report_path = temp / 'combined.tsv'
        if not jobs:
            report_path.write_text('segment\tstrand\treference_gene_index\trank\tgenes_in_path\tgene_in_path\t'
                                   'reference_start_1based\treference_end_1based\tcandidate_start_1based\t'
                                   'candidate_end_1based\tfixed_interval_score\treference_fixed_interval_score\t'
                                   'score_delta\tactual_start_score\tintron_chain\treference_intron_chain\t'
                                   'identical_to_reference\tclassification\twindow_start_1based\t'
                                   'window_end_1based\tinitial_inter_len\n')
        for index, job in enumerate(jobs):
            begin = time.monotonic()
            coords, segment = coordinates(job)
            work = temp / f'decode-{index}'
            work.mkdir()
            scores = {key: job[key] for key in ('emissions', 'gt', 'ag', 'atg', 'stop')}
            raw = work/'paths.gff'
            report = raw.with_suffix('.tsv')
            command = [core, job['fasta'], *(scores[key] for key in ('emissions','gt','ag','atg','stop')),
                       '--k', args.k, '--flank', args.flank, '--threads', threads, '--output', raw, '--report', report]
            if job['log']:
                command.extend(['--viterbi-log', job['log']])
            run_process(command, work/'decoder.log', work)
            converted = work/'converted.gff3'
            stats = convert(raw, segment, job['strand'], converted, work/'stats.json', **coords,
                            k=args.k, short_cds_filter=args.short_cds_filter,
                            overlap=args.segment_overlap, margin=args.segment_margin, best_gff=job['best_gff'])
            with report.open() as handle, report_path.open('w' if index == 0 else 'a') as report_out:
                header = handle.readline()
                if index == 0:
                    report_out.write('segment\tstrand\t' + header)
                for line in handle:
                    report_out.write(f'{segment}\t{job["strand"]}\t{line}')
            # SQLite bounds memory while sorting complete transcript blocks.
            # Child feature order is retained, including for reverse coordinates.
            block, key = [], None
            with converted.open() as handle:
                for line in handle:
                    if line.startswith('#'):
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if fields[2] in ('transcript', 'mRNA'):
                        if block:
                            database.execute('INSERT INTO models VALUES (?,?,?,?,?,?)', (*key, ''.join(block)))
                        key = (fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8].split(';')[0])
                        block = []
                    block.append(line)
                if block:
                    database.execute('INSERT INTO models VALUES (?,?,?,?,?,?)', (*key, ''.join(block)))
            database.commit()
            source = {k:v for k,v in job.items() if k != 'source_inputs'}
            if job.get('source_inputs'):
                source = {k:v for k,v in source.items() if k not in FILE_KEYS}
            source.update(job.get('source_inputs', {}))
            metadata['jobs'].append(dict(inputs={k:str(v) for k,v in source.items() if v is not None},
                                          coordinates=coords, statistics=stats,
                                          seconds=time.monotonic()-begin))
        staged = temp/'result.gff3'
        with staged.open('w') as handle:
            handle.write('##gff-version 3\n')
            for (body,) in database.execute('SELECT body FROM models ORDER BY seqid,start,finish,strand,id'):
                handle.write(body)
        database.close()
        metadata['seconds'] = time.monotonic()-started
        meta_path = temp/'metadata.json'
        meta_path.write_text(json.dumps(metadata, indent=2)+'\n')
        # Each file is replaced atomically, with the main GFF published last.
        os.replace(report_path, targets[1])
        os.replace(meta_path, targets[2])
        os.replace(staged, targets[0])
    print(f'Wrote {output}', file=sys.stderr)


def main():
    args = parser().parse_args()
    try:
        decode(args)
    except (ValueError, OSError, KeyError, ArithmeticError, csv.Error) as error:
        print(f'altann: error: {error}', file=sys.stderr)
        return 2
    return 0
