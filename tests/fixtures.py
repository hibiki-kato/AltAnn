"""Small, deterministic UniAnn inputs with two independent splice choices."""

from pathlib import Path
import csv


def write_fixture(directory: Path, splice_shift: int = 0) -> dict[str, Path]:
    """Write oriented scores and sequence; all positions in inputs are zero based.

    Each gene has a strong coding region and two acceptors three bases apart.
    Both introns preserve frame and exceed the model's 40-base minimum.
    The loci are separated so a 50-base flank cannot combine them.
    """
    directory.mkdir(parents=True, exist_ok=True)
    length = 1200
    sequence = list("C" * length)
    motifs = {
        "atg": ([60, 660], "ATG"),
        "stop": ([450, 1050], "TAA"),
        "gt": ([180 + splice_shift, 780 + splice_shift], "GT"),
        "ag": ([p + splice_shift for p in (238, 241, 838, 841)], "AG"),
    }
    paths = {name: directory / f"{name}.txt" for name in motifs}
    for name, (positions, motif) in motifs.items():
        for position in positions:
            sequence[position:position + len(motif)] = motif
        paths[name].write_text("".join(f"{p}\t2\n" for p in positions))
    paths["fasta"] = directory / "fixture.fa"
    paths["fasta"].write_text(">fixture\n" + "".join(sequence) + "\n")
    paths["emissions"] = directory / "emissions.txt"
    with paths["emissions"].open("w") as stream:
        for position in range(length):
            coding = any(start <= position < start + 393 for start in (60, 660))
            intronic = any(start + splice_shift <= position < start + splice_shift + 60
                           for start in (180, 780))
            values = [-8, -8, -8, -8, -8]
            values[0 if not coding else 1] = 2
            if intronic:
                values[4] = 3
            stream.write(str(position) + "\t" + "\t".join(map(str, values)) + "\n")
    return paths


def ordered_inputs(paths: dict[str, Path]) -> list[str]:
    return [str(paths[key]) for key in ("fasta", "emissions", "gt", "ag", "atg", "stop")]


def write_stranded_fixture(directory: Path, length: int = 1200):
    """One distinct spliced gene on each strand, in genomic score coordinates."""
    directory.mkdir(parents=True, exist_ok=True)
    oriented = [write_fixture(directory / 'legacy-plus'),
                write_fixture(directory / 'legacy-minus', splice_shift=1)]
    sequences = [paths['fasta'].read_text().splitlines()[1] for paths in oriented]
    complement = str.maketrans('ACGT', 'TGCA')
    sequence = sequences[0][:600] + 'C' * (length - 1200) + sequences[1][:600].translate(complement)[::-1]
    paths = {key: directory / ('sequence.fna' if key == 'fasta' else 'out.' +
                              ('ps' if key == 'emissions' else key) + '.txt')
             for key in ('fasta', 'emissions', 'gt', 'ag', 'atg', 'stop')}
    paths['fasta'].write_text('>fixture\n' + sequence + '\n')
    for index, inputs in enumerate(oriented):
        dna = sequence if index == 0 else sequence.translate(complement)[::-1]
        inputs['fasta'].write_text('>fixture\n' + dna + '\n')
        for key in ('gt', 'ag', 'atg', 'stop'):
            rows = inputs[key].read_text().splitlines()
            inputs[key].write_text(''.join(row + '\n' for row in rows if int(row.split()[0]) < 600))
        with inputs['emissions'].open('w') as handle:
            for position in range(length):
                values = [-8, -8, -8, -8, -8]
                values[1 if 60 <= position < 453 else 0] = 2 + index * .125
                if 180 + index <= position < 240 + index:
                    values[4] = 3 + index * .25
                handle.write('\t'.join([str(position), *map(str, values), dna[position]]) + '\n')
    for key in ('emissions', 'gt', 'ag', 'atg', 'stop'):
        with paths[key].open('w') as handle:
            handle.write('# position strand scores\n')
            for strand, inputs in zip(('+', '-'), oriented):
                for line in inputs[key].read_text().splitlines():
                    fields = line.split()
                    position = int(fields[0])
                    if strand == '-':
                        position = length - 1 - position
                    handle.write('\t'.join([str(position), strand, *fields[1:]]) + '\n')
    return paths, *oriented


def write_raw_fixture(directory: Path, length: int = 1400, columns: int = 6):
    """Write raw six-frame PSAURON probabilities and genomic splice-site rows.

    Each strand contains a spliced gene with a coding peak in a different frame.
    The independent inputs use their own oriented FASTA and forward CSV columns.
    """
    directory.mkdir(parents=True, exist_ok=True)
    gene = 'ATG' + 'AAA' * 32 + 'AA' + 'GT' + 'C' * 38 + 'AG'
    gene += 'A' + 'AAA' * 67 + 'TAA'
    sequence = list('C' * length)
    sequence[30:30 + len(gene)] = gene
    minus_genomic_start = length - 650 - len(gene)
    complement = str.maketrans('ACGT', 'TGCA')
    sequence[minus_genomic_start:minus_genomic_start + len(gene)] = gene.translate(complement)[::-1]
    sequence = ''.join(sequence)
    arrays = []
    for frame in range(6):
        values = []
        for index in range(length // 3 + 1):
            coding = ((frame == 0 and 9 <= index <= 130)
                      or (frame == 5 and 216 <= index <= 340))
            values.append(str(.8 + (index % 4) * .06) if coding else '0')
        arrays.append(';'.join(values))
    csv_header = ['sequence', *['summary_' + str(i) for i in range(8)],
                  'forward_all_prob_0', 'forward_all_prob_1', 'forward_all_prob_2',
                  'reverse_all_prob_0', 'reverse_all_prob_1', 'reverse_all_prob_2']
    site_header = ['chrom', 'pos', 'strand', 'type', 'motif', 'prob']
    if columns == 7:
        site_header.insert(-1, 'raw_prob')
    paths = {key: directory / filename for key, filename in
             (('fasta', 'sequence.fna'), ('psauron', 'psauron_score.csv'), ('scores', 'sites.tsv'))}

    def write_inputs(inputs, dna, probabilities, rows):
        inputs['fasta'].write_text('>fixture\n' + dna + '\n')
        with inputs['psauron'].open('w', newline='') as handle:
            handle.write('# PSAURON\n# aggregate\n# probability\n')
            writer = csv.writer(handle)
            writer.writerow(csv_header)
            writer.writerow(['fixture', *['0'] * 8, *probabilities])
        inputs['scores'].write_text('\t'.join(site_header) + '\n' +
                                   ''.join('\t'.join(row) + '\n' for row in rows))

    rows = []
    oriented = []
    for strand, start in (('+', 30), ('-', 650)):
        strand_rows = []
        for position, kind, motif, probability in (
                (start, 'start', 'ATG', '.5'),
                (start + len(gene) - 3, 'stop', 'TAA', '.5'),
                (start + 101, 'donor', 'GT', '.99'),
                (start + 141, 'acceptor', 'AG', '.99')):
            coordinate = position + 1 if strand == '+' else length - position
            tail = ['.01', probability] if columns == 7 else [probability]
            rows.append(['fixture', str(coordinate), strand, kind, motif, *tail])
            strand_rows.append(['fixture', str(position + 1), '+', kind, motif, *tail])
        subdir = directory / ('forward' if strand == '+' else 'reverse')
        subdir.mkdir()
        inputs = {key: subdir / path.name for key, path in paths.items()}
        dna = sequence if strand == '+' else sequence.translate(complement)[::-1]
        probabilities = arrays if strand == '+' else arrays[3:] + arrays[:3]
        write_inputs(inputs, dna, probabilities, strand_rows)
        oriented.append(inputs)
    write_inputs(paths, sequence, arrays, list(reversed(rows)))
    return paths, *oriented
