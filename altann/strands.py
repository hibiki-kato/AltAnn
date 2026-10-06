"""Prepare both genomic strands for the existing oriented native decoder."""
import math


COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                          'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')
DNA = set('ACGTRYSWKMBDHVNacgtryswkmbdhvn')
SCORES = ('emissions', 'gt', 'ag', 'atg', 'stop')


def read_sequence(path):
    header, chunks = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    raise ValueError('Both-strand mode requires exactly one FASTA record')
                header = line.rstrip('\r\n')
                if not header[1:].split():
                    raise ValueError('FASTA header has no sequence identifier')
            else:
                sequence = ''.join(line.split())
                if sequence and header is None:
                    raise ValueError('Sequence found before FASTA header')
                if set(sequence) - DNA:
                    raise ValueError('Invalid DNA symbol in FASTA')
                chunks.append(sequence)
    sequence = ''.join(chunks)
    if header is None or not sequence:
        raise ValueError('Empty FASTA sequence')
    return header, sequence


def split_scores(path, destinations, length, emissions=False):
    """Keep numeric text exact; minus frame labels already refer to RC frames."""
    seen = {strand: bytearray(length) for strand in ('+', '-')}
    counts = {'+': 0, '-': 0}
    with path.open() as source, destinations['+'].open('w') as plus, destinations['-'].open('w') as minus:
        handles = {'+': plus, '-': minus}
        for number, line in enumerate(source, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split()
            expected = (7, 8) if emissions else (3,)
            if len(fields) not in expected or fields[1] not in handles:
                raise ValueError(f'Expected position, strand, and {"five emission scores" if emissions else "transition score"} in {path}:{number}')
            try:
                position = int(fields[0])
            except ValueError as error:
                raise ValueError(f'Invalid score position in {path}:{number}') from error
            strand = fields[1]
            if not 0 <= position < length or seen[strand][position]:
                raise ValueError(f'Invalid or duplicate {strand} score position in {path}:{number}')
            values = fields[2:7] if emissions else fields[2:3]
            try:
                finite = all(math.isfinite(float(value)) for value in values)
            except ValueError:
                finite = False
            if not finite:
                raise ValueError(f'Expected finite scores in {path}:{number}')
            if emissions and len(fields) == 8 and (len(fields[7]) != 1 or fields[7] not in DNA):
                raise ValueError(f'Invalid optional emission nucleotide in {path}:{number}')
            seen[strand][position] = 1
            counts[strand] += 1
            oriented = position if strand == '+' else length - position - 1
            handles[strand].write('\t'.join([str(oriented), *fields[2:]]) + '\n')
    if emissions:
        for strand in ('+', '-'):
            if counts[strand] != length:
                raise ValueError(f'Missing {strand} emission positions in {path}: expected {length}, found {counts[strand]}')


def split_best(path, destinations):
    with path.open() as source, destinations['+'].open('w') as plus, destinations['-'].open('w') as minus:
        handles = {'+': plus, '-': minus}
        for line in source:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9 or fields[6] not in handles:
                raise ValueError(f'Combined best GFF requires nine columns and + or - strands: {path}')
            handles[fields[6]].write(line.rstrip('\r\n') + '\n')


def prepare_both(job, directory):
    directory.mkdir()
    header, sequence = read_sequence(job['fasta'])
    length = len(sequence)
    directories = {strand: directory / strand for strand in ('plus', 'minus')}
    jobs = []
    for strand, oriented in (('plus', sequence), ('minus', sequence.translate(COMPLEMENT)[::-1])):
        work = directories[strand]
        work.mkdir()
        fasta = work / 'sequence.fna'
        fasta.write_text(header + '\n' + oriented + '\n')
        prepared = dict(job, fasta=fasta, strand=strand, log=None,
                        source_inputs={key: job[key] for key in ('fasta', *SCORES, 'best_gff') if job.get(key)})
        for key in SCORES:
            prepared[key] = work / f'{key}.txt'
        prepared['best_gff'] = work / 'best.gff3' if job.get('best_gff') else None
        jobs.append(prepared)
    for key in SCORES:
        split_scores(job[key], {'+': jobs[0][key], '-': jobs[1][key]}, length, emissions=key == 'emissions')
    if job.get('best_gff'):
        split_best(job['best_gff'], {'+': jobs[0]['best_gff'], '-': jobs[1]['best_gff']})
    return jobs
