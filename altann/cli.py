"""Read existing UniAnn inputs, invoke the native decoder, and export GFF3.

The native decoder owns all dynamic programming. This module only prepares
files, maps coordinates, and records provenance. A UniAnn log supplies scaling
parameters; its rounded DP trace is never used as decoder input.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time

from . import __version__
from .convert import SEGMENT, convert

# Per-frame chromosome scores routinely exceed csv's 128 KiB default.
# Some platforms expose a smaller C long than Python's integer range.
_csv_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_csv_limit)
        break
    except OverflowError:
        _csv_limit //= 10

ROOT = Path(__file__).resolve().parent.parent
FILE_KEYS = ('fasta', 'emissions', 'gt', 'ag', 'atg', 'stop', 'psauron', 'sites', 'log', 'best_gff')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=f'AltAnn {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    d = sub.add_parser('decode', help='Decode existing UniAnn score files or original inputs.')
    d.add_argument('--input', type=Path, help='Directory containing exactly one sequence and its scores.')
    d.add_argument('--manifest', type=Path, help='TSV inventory; relative paths are resolved beside the manifest.')
    for key in FILE_KEYS:
        d.add_argument('--' + key.replace('_', '-'), type=Path)
    d.add_argument('--strand', choices=('plus', 'minus'), help='Orientation of the supplied FASTA; minus must already be reverse complemented.')
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


def jobs_from_args(args):
    if args.manifest:
        if args.input or any(getattr(args, key) for key in FILE_KEYS):
            raise ValueError('--manifest cannot be combined with --input or individual file arguments')
        with args.manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if not rows:
            raise ValueError('Manifest contains no jobs')
        base = args.manifest.resolve().parent
    else:
        rows = [{key: getattr(args, key) for key in FILE_KEYS}]
        rows[0].update(strand=args.strand, seqid=args.seqid, offset=args.offset,
                       sequence_length=args.sequence_length, score_dir=args.input)
        base = Path.cwd()
    for row in rows:
        for key in FILE_KEYS + ('score_dir',):
            value = row.get(key)
            row[key] = (base / str(value)).resolve() if value else None
        directory = row['score_dir'] or (row['fasta'].parent if row['fasta'] else None)
        if directory is None:
            raise ValueError('Supply --input, --fasta, or a manifest containing fasta paths')
        if not directory.is_dir():
            raise ValueError(f'Input directory does not exist: {directory}')
        if not row['fasta']:
            row['fasta'] = unique(directory, ('*.fa', '*.fasta', '*.fna'), 'FASTA', True)
        for key, filename in [('emissions','out.ps.txt'), ('gt','out.gt.txt'), ('ag','out.ag.txt'), ('atg','out.atg.txt'), ('stop','out.stop.txt')]:
            if not row[key] and (directory / filename).is_file():
                row[key] = directory / filename
        if not row['emissions'] and not row['psauron']:
            row['psauron'] = unique(directory, ('psauron_score.csv',), 'PSAURON scores', True)
        if not all(row[k] for k in ('gt','ag','atg','stop')) and not row['sites']:
            row['sites'] = unique(directory, ('sites*.tsv',), 'site scores', True)
        if row['sites'] and not row['log']:
            row['log'] = unique(directory, ('uniann.log', '*.uniann.log', 'run.log'), 'UniAnn log')
        if not row['best_gff']:
            row['best_gff'] = unique(directory, ('*.uniann.gff',), 'UniAnn best GFF')
        for key in FILE_KEYS:
            if row[key] and not row[key].is_file():
                raise ValueError(f'Missing {key} file: {row[key]}')
        if not row.get('strand'):
            # Only explicit conventional filenames or directory tokens imply orientation.
            tokens = re.split(r'[._-]', directory.name.lower())
            if row['fasta'].name.lower().endswith('.rc.fa') or 'minus' in tokens:
                row['strand'] = 'minus'
            elif 'plus' in tokens:
                row['strand'] = 'plus'
            else:
                raise ValueError('Input orientation is ambiguous; specify --strand plus or minus')
        if row['strand'] not in ('plus', 'minus'):
            raise ValueError(f"Invalid strand: {row['strand']}")
    return rows


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
                length += len(line.strip())
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


def prepare_sites(job, directory):
    """Reproduce UniAnn's plus-oriented site score transform and serialization."""
    # The shell wrapper obtains exp(1) through Perl's 15-digit string output.
    multiplier, factor = float(format(math.exp(1), '.15g')), None
    if job['log']:
        pattern = re.compile(r'Multiplier is ([^,\s]+), Factor is (\S+)')
        values = set()
        with job['log'].open() as handle:
            for line in handle:
                match = pattern.search(line)
                if match:
                    values.add(tuple(map(float, match.groups())))
        if len(values) != 1:
            raise ValueError('Log must contain one unambiguous Multiplier/Factor pair')
        multiplier, factor = values.pop()
    maximum = 0.0
    with job['fasta'].open() as handle:
        sequence_id = next(line[1:].split()[0] for line in handle if line.startswith('>'))
    # UniAnn uses column seven for factor estimation, but the last column for
    # final site scores. Preserve this distinction for files with extra fields.
    def rows():
        with job['sites'].open() as handle:
            for line in handle:
                fields = line.split()
                if not fields or line.startswith('#'):
                    continue
                if len(fields) < 6:
                    raise ValueError('Site table requires at least six whitespace-separated columns')
                # The saved UniAnn site table normally starts with this header.
                if fields[:4] == ['chrom', 'pos', 'strand', 'type']:
                    continue
                if fields[0] != sequence_id:
                    raise ValueError(f'Site sequence identifier {fields[0]} does not match FASTA {sequence_id}')
                if fields[2] == '+':
                    yield fields
    if factor is None:
        for fields in rows():
            if fields[3] == 'donor':
                value = float(fields[6] if len(fields) > 6 else fields[5])
                maximum = max(maximum, math.log(value * multiplier + 1e-10))
        if maximum <= 0:
            raise ValueError('Cannot infer scaling: no donor has a positive transformed score')
        factor = int(1000 / maximum + .5)
    if not math.isfinite(multiplier) or not math.isfinite(factor) or multiplier <= 0 or factor <= 0:
        raise ValueError('Multiplier and factor must be finite and positive')
    names = {'donor':'gt', 'acceptor':'ag', 'start':'atg', 'stop':'stop'}
    handles = {key: (directory / f'out.{key}.txt').open('w') for key in names.values()}
    try:
        for fields in rows():
            if fields[3] not in names:
                continue
            probability = float(fields[-1])
            if not math.isfinite(probability) or not 0 <= probability <= 1:
                raise ValueError('Site probabilities must be finite and between zero and one')
            score = math.log(probability * multiplier + 1e-10) * factor
            handles[names[fields[3]]].write(f'{int(fields[1])-1}\t{score:.15g}\n')
    finally:
        for handle in handles.values():
            handle.close()
    return {key: directory / f'out.{key}.txt' for key in names.values()}, dict(multiplier=multiplier, factor=factor)


def validate_psauron(path, fasta):
    """Reject missing or malformed PSAURON frames before legacy preprocessing.

    UniAnn's preprocessing format has four preamble lines and uses columns
    10--12 for the three forward frames. Reverse inputs are already oriented.
    """
    with fasta.open() as handle:
        sequence_id = next(line[1:].split()[0] for line in handle if line.startswith('>'))
    found = 0
    with path.open() as handle:
        for _ in range(4):
            if not handle.readline():
                raise ValueError('PSAURON input is missing its four-line preamble')
        for row in csv.reader(handle):
            if not row or row[0] != sequence_id:
                continue
            found += 1
            if len(row) < 15:
                raise ValueError('PSAURON input requires all frame columns (run PSAURON with -a)')
            for frame in row[9:12]:
                if not frame:
                    raise ValueError('Empty PSAURON forward frame')
                for value in frame.rstrip(';').split(';'):
                    value = float(value)
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        raise ValueError('PSAURON probabilities must be finite and between zero and one')
    if found != 1:
        raise ValueError(f'Expected one PSAURON row for FASTA identifier {sequence_id}; found {found}')


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
    jobs = jobs_from_args(args)
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
    inputs = {Path(v).resolve() for job in jobs for k,v in job.items() if k in FILE_KEYS and v}
    if args.manifest:
        inputs.add(args.manifest.resolve())
    if any(t in inputs or (t.exists() and any(t.samefile(p) for p in inputs)) for t in targets):
        raise ValueError('Output paths must not overwrite input files')
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(version=__version__, command=sys.argv, threads=threads, jobs=[])
    with tempfile.TemporaryDirectory(prefix='.altann-', dir=output.parent) as temp:
        temp = Path(temp)
        database = sqlite3.connect(temp / 'merge.sqlite')
        database.execute('CREATE TABLE models (seqid TEXT, start INTEGER, finish INTEGER, strand TEXT, id TEXT, body TEXT)')
        report_path = temp / 'combined.tsv'
        report_header = None
        identities = set()
        with report_path.open('w') as report_out:
            for index, job in enumerate(jobs):
                begin = time.monotonic()
                coords, segment = coordinates(job)
                identity = (segment, job['strand'])
                if identity in identities:
                    raise ValueError(f'Duplicate segment/strand identifier: {identity}')
                identities.add(identity)
                work = temp / str(index)
                work.mkdir()
                scores = {key: job[key] for key in ('emissions','gt','ag','atg','stop')}
                scaling = None
                if not scores['emissions']:
                    validate_psauron(job['psauron'], job['fasta'])
                    run_process(['perl', ROOT/'vendor'/'preprocess_psauron_scores.pl', job['fasta'], job['psauron']], work/'preprocess.log', work)
                    scores['emissions'] = work/'out.ps.txt'
                if not all(scores[key] for key in ('gt','ag','atg','stop')):
                    generated, scaling = prepare_sites(job, work)
                    for key, path in generated.items():
                        scores[key] = scores[key] or path
                raw = work/'paths.gff'
                report = raw.with_suffix('.tsv')
                command = [core, job['fasta'], *(scores[key] for key in ('emissions','gt','ag','atg','stop')),
                           '--k', args.k, '--flank', args.flank, '--threads', threads, '--output', raw, '--report', report]
                run_process(command, work/'decoder.log', work)
                converted = work/'converted.gff3'
                stats = convert(raw, segment, job['strand'], converted, work/'stats.json', **coords,
                                k=args.k, short_cds_filter=args.short_cds_filter,
                                overlap=args.segment_overlap, margin=args.segment_margin, best_gff=job['best_gff'])
                with report.open() as handle:
                    header = handle.readline().rstrip('\n')
                    if report_header is None:
                        report_header = header
                        report_out.write('segment\tstrand\t'+header+'\n')
                    elif header != report_header:
                        raise ValueError('Inconsistent native report columns')
                    for line in handle:
                        report_out.write(f'{segment}\t{job["strand"]}\t{line}')
                # SQLite bounds memory while sorting complete transcript blocks
                # across chromosomes and segments. Child feature order is retained.
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
                metadata['jobs'].append(dict(inputs={k:str(v) for k,v in job.items() if v is not None},
                                              coordinates=coords, statistics=stats, scaling=scaling,
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
