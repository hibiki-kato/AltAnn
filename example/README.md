# S. pombe tutorial

Run AltAnn on chromosome III of *Schizosaccharomyces pombe*
(`NC_003421.2`) using the supplied FASTA and emission and transition scores.

## Before you start

Install AltAnn 0.2.0 or later and Python 3.9 or later using the
[installation instructions](../README.md#install-a-release).
No additional Python packages, scoring software, or GPU are needed.

Run the commands below from the repository root. If you installed a release
package, replace `bin/altann` with the extracted package's launcher path.

The inputs are ready to use:

| Directory | FASTA | Scores |
| --- | --- | --- |
| `example/data/NC_003421.2/plus/` | `sequence.fna`: original sequence | Scores for the original sequence |
| `example/data/NC_003421.2/minus/` | `sequence.rc.fna`: reverse complement | Scores for the reverse complement |

Each directory contains the emission file `out.ps.txt` and transition files
`out.gt.txt`, `out.ag.txt`, `out.atg.txt`, and `out.stop.txt`.

## 1. Decode the forward input

```sh
bin/altann decode --input example/data/NC_003421.2/plus \
  --rerun-viterbi --k 10 --threads 4 --output example/results/chr3.plus.gff3
```

`--rerun-viterbi` computes the best path from the supplied scores before
searching for alternatives. This example uses that mode because a baseline
GFF and Viterbi log are not included.

`--k 10` retains up to ten paths per state during the search. Use `--k 5` to
retain fewer paths. K is not the number of exported isoforms.
`--threads 4` uses four threads; use `--threads auto` to use available CPUs.

## 2. Decode the reverse input

**`--reverse` requires an already reverse-complemented FASTA and matching
scores.** It does not transform the input FASTA for you.

The supplied `minus/sequence.rc.fna` is already the reverse complement of
`plus/sequence.fna`, and the scores in `minus/` match it. Use these files
as supplied, without converting them again.

```sh
bin/altann decode --input example/data/NC_003421.2/minus \
  --reverse --rerun-viterbi --k 10 --threads 4 \
  --output example/results/chr3.minus.gff3
```

AltAnn maps the reverse output back to the original chromosome coordinates
and uses strand `-`. The forward output uses strand `+`.
Each invocation processes one direction.

For your own reverse input, first reverse complement the FASTA, then generate
scores for that sequence. Any supplied UniAnn GFF and Viterbi log must also
match that orientation. See [input formats](../docs/format.md) for details.

## 3. Read the results

The results are written to `example/results/`:

| GFF3 file | Expected transcripts with K=10 |
| --- | --- |
| `chr3.plus.gff3` | 1,114: 385 references and 729 alternatives |
| `chr3.minus.gff3` | 1,108: 402 references and 706 alternatives |

Each run also writes a `.gff3.tsv` diagnostics table and a `.gff3.json` record
of the run settings and inputs.

The GFF3 contains transcripts, exons, and CDS. Column 6 is `.`; scores are
stored in transcript attributes:

- `kbest_score`: path score.
- `kbest_reference_score`: reference path score.
- `kbest_delta`: score difference from the reference.
- `kbest_rank`: path rank; rank 0 identifies the reference.
- `kbest_origin`: reference or alternative model.

Path scores are comparable within a locus and are not probabilities.
Both GFF3 files use the original chromosome coordinates, so they can be
combined with a downstream GFF tool when a single annotation file is needed.
