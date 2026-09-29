#!/usr/bin/env python3
"""Prepare strand-oriented S. pombe inputs using only the Python standard library."""

import argparse
import gzip
import shutil
from pathlib import Path


CHROMOSOMES = ("NC_003424.3", "NC_003423.3", "NC_003421.2")
STRANDS = ("plus", "minus")
SCORE_FILES = ("psauron_score.csv", "sites.tsv")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def read_chromosomes(archive):
    """Read nuclear sequences; the mitochondrial record is not used by this example."""
    sequences = {}
    identifier = None
    with gzip.open(archive, "rt", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                if identifier in CHROMOSOMES:
                    if identifier in sequences:
                        raise ValueError(f"Duplicate FASTA record: {identifier}")
                    sequences[identifier] = []
            elif identifier in CHROMOSOMES:
                sequences[identifier].append(line.strip().upper())
    missing = [name for name in CHROMOSOMES if name not in sequences]
    if missing:
        raise ValueError("Missing nuclear chromosomes: " + ", ".join(missing))
    for name, chunks in sequences.items():
        sequence = "".join(chunks)
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Empty sequence or unsupported DNA symbols in {name}")
        sequences[name] = sequence
    return sequences


def write_fasta(path, identifier, sequence):
    """Keep the original accession on both strands for genomic GFF coordinates."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(f">{identifier}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset:offset + 80] + "\n")


def prepare(data_dir, work_dir, genome_only=False):
    genome = data_dir / "genome.fna.gz"
    required = [genome]
    if not genome_only:
        required.extend(
            data_dir / name / strand / (filename + ".gz")
            for name in CHROMOSOMES
            for strand in STRANDS
            for filename in SCORE_FILES
        )
    # Check the entire input set before creating any jobs. In particular, a
    # genome alone must never be presented as a ready-to-decode example.
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError("Missing input archives:\n  " + "\n  ".join(missing))

    sequences = read_chromosomes(genome)
    rows = ["fasta\tstrand\tscore_dir\n"]
    for name in CHROMOSOMES:
        for strand in STRANDS:
            job = work_dir / name / strand
            sequence = sequences[name]
            filename = "sequence.fa"
            if strand == "minus":
                sequence = sequence.translate(COMPLEMENT)[::-1]
                filename = "sequence.rc.fa"
            write_fasta(job / filename, name, sequence)
            if genome_only:
                continue
            for score_file in SCORE_FILES:
                archive = data_dir / name / strand / (score_file + ".gz")
                # Score tables can be much larger than the genome. Copy them
                # in bounded chunks rather than reading them into memory.
                with gzip.open(archive, "rb") as source:
                    with (job / score_file).open("wb") as destination:
                        shutil.copyfileobj(source, destination, length=1024 * 1024)
            relative_job = job.relative_to(work_dir).as_posix()
            rows.append(f"{relative_job}/{filename}\t{strand}\t{relative_job}\n")
    if not genome_only:
        (work_dir / "inputs.tsv").write_text("".join(rows), encoding="ascii")


def main():
    example = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=example / "data")
    parser.add_argument("--work-dir", type=Path, default=example / "work")
    parser.add_argument(
        "--genome-only", action="store_true",
        help="prepare FASTA files for score regeneration; do not write inputs.tsv",
    )
    args = parser.parse_args()
    try:
        prepare(args.data_dir, args.work_dir, args.genome_only)
    except (OSError, EOFError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")
    label = "FASTA files" if args.genome_only else "FASTA files, scores, and inputs.tsv"
    print(f"Prepared {label} in {args.work_dir}")


if __name__ == "__main__":
    main()
