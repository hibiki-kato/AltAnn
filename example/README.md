# S. pombe example

Practice AltAnn on chromosome III of *Schizosaccharomyces pombe*
(`NC_003421.2`, 2,452,883 bases, RefSeq assembly `GCF_000002945.2` / ASM294v3).
The forward sequence and its reverse complement each have their own prepared
emission and transition score tables. No scoring software or GPU is needed.

The compressed inputs total 34.7 MB. Processed scores for all three chromosomes
would exceed the 100 MB example budget, so this tutorial uses chromosome III.
Checksums are in `SHA256SUMS`; data generation and model details are in
`provenance.json`. Those scoring tools were used to prepare the files and are
not AltAnn dependencies.

## Prepare the input files

Use AltAnn 0.2.0 or later and Python 3.9 or later. Follow the repository's
[build instructions](../README.md#build-from-source-developers), or download a
[release package](https://github.com/hibiki-kato/AltAnn/releases).
Run these commands from the repository root:

```sh
python3 example/prepare.py
```

This only decompresses the supplied files into `example/work/inputs/`.
It does not calculate scores or reverse sequences. Each orientation contains
one FASTA and `out.ps.txt`, `out.gt.txt`, `out.ag.txt`, `out.atg.txt`, and
`out.stop.txt`. The latter four files are position-specific transition scores.
The legacy emission format has N, E0, E1, E2, and a shared I score for the three
intron states, giving the seven-state UniAnn model.

## Decode the forward input

```sh
bin/altann decode --input example/work/inputs/NC_003421.2/plus \
  --rerun-viterbi --k 10 --threads 4 --output example/results/chr3.plus.gff3
```

The example omits the original large DP/BT logs and GFF. `--rerun-viterbi`
explicitly reconstructs the best path from the prepared scores before K-best
search. Normal operation can reuse a matching UniAnn GFF and log instead;
see [input formats](../docs/format.md).

## Decode the reverse input separately

```sh
bin/altann decode --input example/work/inputs/NC_003421.2/minus \
  --reverse --rerun-viterbi --k 10 --threads 4 \
  --output example/results/chr3.minus.gff3
```

`--reverse` tells AltAnn that this sequence has already been reverse
complemented. All scores use that sequence's coordinates. AltAnn leaves these
inputs unchanged and maps only the output GFF to the original chromosome,
using the negative strand. One run always processes one direction.

For release packages, replace `bin/altann` with the extracted package's launcher
path. `--threads auto` uses available CPUs; `--k 5` retains fewer local paths.
K is not a count of exported isoforms: duplicate intron chains are removed,
and a path can contain more than one gene.

## Read the output

With K=10, the forward run exports 1,114 transcripts (385 references and 729
alternatives), and the reverse run exports 1,108 (402 references and 706
alternatives). Both outputs were verified against a separate run importing the
original UniAnn DP/BT logs; see [validation](../docs/validation.md).

Each run writes its GFF, a `.gff3.tsv` diagnostics table, and `.gff3.json`
provenance file. `example/work/` and `example/results/` are ignored by Git.
The GFF contains transcripts, exons, and CDS. Column 6 is `.`; transcript
attributes carry `kbest_score`, `kbest_reference_score`, `kbest_delta`,
`kbest_rank`, and `kbest_origin`. Rank 0 is the UniAnn reference. Path scores
are comparable within a locus and are not probabilities.

Both GFF files use the original chromosome coordinates and different strand
signs. They can be combined by a downstream GFF tool if a single annotation
file is needed. AltAnn does not decode both directions in the same invocation.
