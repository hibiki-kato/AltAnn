#!/usr/bin/env python3
"""Stage the bundled original FASTA and both-strand probability tables."""
import gzip
import os
from pathlib import Path
import shutil
import tempfile


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'data' / 'NC_003421.2'
OUTPUT = ROOT / 'work' / 'NC_003421.2'


def prepare():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # Stage complete files before replacing the reusable example inputs.
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=OUTPUT) as temp:
        temp = Path(temp)
        shutil.copyfile(SOURCE / 'plus' / 'sequence.fna', temp / 'sequence.fna')
        for name in ('psauron_score.csv', 'sites.tsv'):
            with gzip.open(SOURCE / 'raw' / (name + '.gz'), 'rb') as source:
                with (temp / name).open('wb') as destination:
                    shutil.copyfileobj(source, destination)
        for path in temp.iterdir():
            os.replace(path, OUTPUT / path.name)
    print(f'Prepared both-strand inputs: {OUTPUT}')


if __name__ == '__main__':
    prepare()
