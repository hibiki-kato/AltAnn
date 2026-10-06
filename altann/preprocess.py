"""Convert original UniAnn probability inputs into oriented decoder scores."""
from contextlib import ExitStack
import csv
import math
import re
import sys

from .strands import COMPLEMENT, SCORES, read_sequence, split_best


STOP = -1e6
SITE_TYPES = {'donor': 'gt', 'acceptor': 'ag', 'start': 'atg', 'stop': 'stop'}
DEFAULT_MULT = float('2.71828182845905')


def probability(value, label):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'Invalid probability in {label}: {value!r}') from error
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        raise ValueError(f'Probabilities must be finite and between zero and one in {label}')
    return parsed


def read_frames(path, seqid, both):
    """Read PSAURON's original forward and reverse arrays without reordering."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            limit //= 10
    frames = None
    with path.open(newline='') as handle:
        for _ in range(4):
            if not handle.readline():
                raise ValueError('PSAURON input is missing its four-line preamble')
        try:
            for row in csv.reader(handle):
                if not row or row[0] != seqid:
                    continue
                if frames is not None:
                    raise ValueError(f'Duplicate PSAURON scores for {seqid}')
                if len(row) < (15 if both else 12):
                    raise ValueError('Missing PSAURON frame columns; run PSAURON with -a for both strands')
                frames = []
                for column in range(9, 15 if both else 12):
                    values = row[column].split(';')
                    if column < 12:
                        while values and values[-1] == '':
                            values.pop()
                    if not values or any(value == '' for value in values):
                        raise ValueError(f'Empty PSAURON frame {column - 9} for {seqid}')
                    frames.append([probability(value, f'{path}, frame {column - 9}') for value in values])
        except csv.Error as error:
            raise ValueError(f'Malformed PSAURON CSV: {path}: {error}') from error
    if frames is None:
        raise ValueError(f'Missing PSAURON scores for {seqid}')
    return {'+': frames[:3], '-': frames[3:]}


def write_emissions(sequence, frames, path):
    """Preserve UniAnn's stop insertion, carried scores, and two-decimal output."""
    stops = {match.start() for match in re.finditer('TAA|TAG|TGA', sequence.upper())}
    insertions = [{}, {}, {}]
    counts = [0, 0, 0]
    # UniAnn's elsif chain records only the first stopped frame in each group.
    for index, position in enumerate(range(0, len(sequence) - 3, 3)):
        for frame in range(3):
            if position + frame in stops:
                before = index - counts[frame]
                insertions[frame][before] = insertions[frame].get(before, 0) + 1
                counts[frame] += 1
                break
    expanded = []
    for frame, values in enumerate(frames):
        scores = []
        # The upstream insertion loop omits markers beyond the supplied array.
        for index, value in enumerate(values):
            scores.extend([STOP] * insertions[frame].get(index, 0))
            scores.append(math.log(value * 50 + 1e-6) / math.log(50))
        expanded.append(scores)
    previous = [0.0, 0.0, 0.0]
    with path.open('w') as handle:
        for position, base in enumerate(sequence):
            if base not in 'ACGTacgt':
                scores = [0, STOP, STOP, STOP, 0]
            else:
                coding = []
                for frame in range(3):
                    index = (position - frame) // 3
                    value = expanded[frame][index] if position >= frame and index < len(expanded[frame]) else 0.0
                    if value == STOP and (position - frame) % 3 in (0, 1):
                        value = previous[frame]
                    coding.append(value)
                ordered = sorted(coding)
                intergenic = 0.15 - (ordered[2] - ordered[1]) if ordered[2] > 0 else 4
                intron = 1.0 - ordered[2]
                if STOP in coding:
                    intergenic, intron = 8, -2
                scores = [intergenic, *coding, intron]
                for frame, value in enumerate(coding):
                    if value > STOP:
                        previous[frame] = value
            handle.write(f'{position}\t' + '\t'.join(f'{value:.2f}' for value in scores) + f'\t{base}\n')


def site_rows(path, seqid, length, strands):
    """Stream matching original-coordinate rows while validating probabilities."""
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split()
            if fields[:4] == ['chrom', 'pos', 'strand', 'type'] or fields[0] != seqid:
                continue
            label = f'{path}:{number}'
            if len(fields) not in (6, 7):
                raise ValueError(f'Expected six or seven site columns in {label}')
            if fields[2] not in ('+', '-'):
                raise ValueError(f'Invalid site strand in {label}')
            if fields[2] not in strands:
                continue
            try:
                position = int(fields[1])
            except ValueError as error:
                raise ValueError(f'Invalid site position in {label}') from error
            if not 1 <= position <= length:
                raise ValueError(f'Site position outside FASTA sequence in {label}')
            if fields[3] not in SITE_TYPES:
                raise ValueError(f'Unknown site type in {label}: {fields[3]}')
            value = probability(fields[5], label)
            if len(fields) == 7:
                value = probability(fields[6], label)
            yield fields[2], position, SITE_TYPES[fields[3]], value


def prepare_raw(job, directory, both=True):
    """Prepare raw UniAnn inputs and skip strands without a qualifying donor."""
    header, sequence = read_sequence(job['fasta'])
    seqid, length = header[1:].split()[0], len(sequence)
    frames = read_frames(job['psauron'], seqid, both)
    try:
        multiplier = float(job.get('mult', DEFAULT_MULT))
    except (TypeError, ValueError) as error:
        raise ValueError('Multiplier must be between 1 and 100') from error
    if not math.isfinite(multiplier) or not 1 <= multiplier <= 100:
        raise ValueError('Multiplier must be between 1 and 100')
    strands = ('+', '-') if both else ('+',)
    maximum = {strand: 0.0 for strand in strands}
    for strand, _, kind, value in site_rows(job['scores'], seqid, length, strands):
        if kind == 'gt':
            maximum[strand] = max(maximum[strand], math.log(value * multiplier + 1e-10))
    directory.mkdir()
    jobs, prepared = [], {}
    for strand in strands:
        if maximum[strand] <= 0:
            print(f'Skipping {strand} strand: no donor with a positive rescaled score', file=sys.stderr)
            continue
        orientation = ('plus' if strand == '+' else 'minus') if both else job.get('strand', 'plus')
        work = directory / ('plus' if strand == '+' else 'minus')
        work.mkdir()
        oriented = sequence if strand == '+' else sequence.translate(COMPLEMENT)[::-1]
        item = dict(job, fasta=work / 'sequence.fna', strand=orientation, log=None,
                    source_inputs={key: job[key] for key in ('fasta', 'psauron', 'scores', 'best_gff') if job.get(key)})
        item.update({key: work / f'out.{"ps" if key == "emissions" else key}.txt' for key in SCORES})
        item['best_gff'] = work / 'best.gff3' if both and job.get('best_gff') else job.get('best_gff')
        item['fasta'].write_text(header + '\n' + oriented + '\n')
        write_emissions(oriented, frames[strand], item['emissions'])
        prepared[strand] = (item, math.floor(1000 / maximum[strand] + .5))
        jobs.append(item)
    with ExitStack() as stack:
        handles = {(strand, kind): stack.enter_context(item[kind].open('w'))
                   for strand, (item, _) in prepared.items() for kind in SITE_TYPES.values()}
        for strand, position, kind, value in site_rows(job['scores'], seqid, length, strands):
            if strand not in prepared:
                continue
            native_position = position - 1 if strand == '+' else length - position
            score = math.log(value * multiplier + 1e-10) * prepared[strand][1]
            handles[strand, kind].write(f'{native_position}\t{score:.15g}\n')
    if both and job.get('best_gff'):
        destinations = {strand: prepared[strand][0]['best_gff'] if strand in prepared else directory / f'skipped-{strand}.gff3'
                        for strand in ('+', '-')}
        split_best(job['best_gff'], destinations)
        # A skipped strand has an empty baseline, so it cannot match supplied features.
        for strand, path in destinations.items():
            if strand not in prepared and path.stat().st_size:
                raise ValueError(f'Supplied GFF contains features for skipped {strand} strand')
    return jobs
