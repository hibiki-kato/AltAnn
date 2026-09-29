#!/usr/bin/env python3
"""Expand prepared S. pombe inputs without transforming sequences or scores."""

import argparse
import gzip
import shutil
from pathlib import Path


def prepare(data_dir, work_dir):
    """Require both independent input sets before extracting the tutorial files."""
    archives = []
    for strand in ("plus", "minus"):
        fasta = "sequence.fa" if strand == "plus" else "sequence.rc.fa"
        names = [fasta] + [f"out.{kind}.txt" for kind in ("ps", "gt", "ag", "atg", "stop")]
        for name in names:
            relative = Path("NC_003421.2") / strand / name
            archive = data_dir / relative.with_suffix(relative.suffix + ".gz")
            if not archive.is_file():
                raise ValueError(f"Missing prepared input: {archive}")
            archives.append((archive, work_dir / relative))
    for archive, destination in archives:
        destination.parent.mkdir(parents=True, exist_ok=True)
        # Score tables are large; bounded copies avoid loading them into RAM.
        with gzip.open(archive, "rb") as source, destination.open("wb") as output:
            shutil.copyfileobj(source, output, length=1024 * 1024)


def main():
    example = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=example / "data")
    parser.add_argument("--work-dir", type=Path, default=example / "work" / "inputs")
    args = parser.parse_args()
    try:
        prepare(args.data_dir, args.work_dir)
    except (OSError, EOFError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(f"Prepared independent forward and reverse inputs in {args.work_dir}")


if __name__ == "__main__":
    main()
