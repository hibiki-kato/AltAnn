# S. pombe example

Practice AltAnn on chromosome III of *Schizosaccharomyces pombe*
(`NC_003421.2`, 2,452,883 bases, RefSeq assembly `GCF_000002945.2` / ASM294v3).
The forward sequence and its reverse complement each have their own prepared
emission and transition score tables. No scoring software or GPU is needed.

The inputs are plain files tracked by Git, without Git LFS. Each file is below
100 MB; the largest emission table is 89.4 MB. The two directions total
206.8 MB. Chromosomes I and II are omitted because each of their full emission
tables exceeds 100 MB. Checksums are in `SHA256SUMS`; data generation and model
details are in `provenance.json`. The scoring tools listed there prepared the
files and are not AltAnn dependencies.

## Locate the input files

Use AltAnn 0.2.0 or later and Python 3.9 or later. Follow the repository's
[build instructions](../README.md#build-from-source-developers), or download a
[release package](https://github.com/hibiki-kato/AltAnn/releases).
Run the commands below from the repository root. No decompression or
preparation step is required.

The forward FASTA is `example/data/NC_003421.2/plus/sequence.fna`; the already
reverse-complemented FASTA is `example/data/NC_003421.2/minus/sequence.rc.fna`.
Each directory also contains
`out.ps.txt`, `out.gt.txt`, `out.ag.txt`, `out.atg.txt`, and `out.stop.txt`.
The latter four files are position-specific transition scores. The legacy
emission format has N, E0, E1, E2, and a shared I score for the three intron
states, giving the seven-state UniAnn model.

## Decode the forward input

```sh
bin/altann decode --input example/data/NC_003421.2/plus \
  --rerun-viterbi --k 10 --threads 4 --output example/results/chr3.plus.gff3
```

The example omits the original large DP/BT logs and GFF. `--rerun-viterbi`
explicitly reconstructs the best path from the prepared scores before K-best
search. Normal operation can reuse a matching UniAnn GFF and log instead;
see [input formats](../docs/format.md).

## Decode the reverse input separately

**For `--reverse`, the input FASTA must already be reverse complemented.**
This requires reversing the sequence and complementing the bases (A with T,
C with G). The flag does not transform the FASTA for you.

In this example, `minus/sequence.rc.fna` is already the reverse complement of
`plus/sequence.fna`; their sequence identifiers and lengths are the same, but
their base sequences differ. The score files in `minus/` were prepared for that
reverse-complemented sequence. Use them as supplied, without converting them
again. For your own data, prepare the reverse-complemented FASTA first, then
generate matching scores and, if used, a UniAnn GFF and Viterbi log.

```sh
bin/altann decode --input example/data/NC_003421.2/minus \
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
