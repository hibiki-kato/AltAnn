# S. pombe tutorial

Run AltAnn on both strands of chromosome III of *Schizosaccharomyces pombe*
(`NC_003421.2`) using the supplied original FASTA, PSAURON CSV, and site probabilities.
The example includes only chromosome III. Chromosomes I and II are not included.

## Before you start

Build the current source or install AltAnn 0.3.0 or later and Python 3.9 or later
using the [installation instructions](../README.md#install-a-release).
No additional Python packages, scoring software, or GPU are needed.

Run the commands below from the repository root. If you installed a release
package, replace `bin/altann` with the extracted package's launcher path.

## 1. Prepare one input set

Stage the original FASTA and decompress the bundled both-strand probability tables:

```sh
python3 example/prepare_stranded_inputs.py
```

This writes `sequence.fna`, `psauron_score.csv`, and `sites.tsv` into
`example/work/NC_003421.2/`. The CSV contains all six frame arrays; the site
probability table uses one-based original coordinates and `+`/`-` rows.
The tables combine the previously scored forward and reverse-complemented
sequences. Preparation does not run PSAURON or generate new probabilities.

## 2. Decode both strands

```sh
bin/altann decode --input example/work/NC_003421.2 -a \
  --k 10 --threads 4 --output example/results/chr3.gff3
```

`-a` (`--all-prob`) uses the original FASTA and the scores for both strands.
AltAnn prepares the reverse complement and emission/transition scores internally,
computes both reference paths, and writes one sorted GFF3 in original coordinates. A baseline
GFF or Viterbi log is not needed for this mode.

`--k 10` retains up to ten paths per state during the search. Use `--k 5` to
retain fewer paths. K is not the number of exported isoforms.
`--threads 4` uses four threads; use `--threads auto` to use available CPUs.

For your own input, supply a PSAURON six-frame CSV and both-strand site
probabilities as described in [input formats](../docs/format.md#both-strand-scores).

## 3. Read the results

`example/results/chr3.gff3` contains 2,222 transcripts with K=10:
787 references and 1,435 alternatives. There are 1,114 plus-strand transcripts
and 1,108 minus-strand transcripts, matching two independent oriented runs.

The run also writes a `.gff3.tsv` diagnostics table containing both strands
and a `.gff3.json` record of the settings and original inputs for each strand.

The GFF3 contains transcripts, exons, and CDS. Column 6 is `.`; scores are
stored in transcript attributes:

- `kbest_score`: path score.
- `kbest_reference_score`: reference path score.
- `kbest_delta`: score difference from the reference.
- `kbest_rank`: path rank; rank 0 identifies the reference.
- `kbest_origin`: reference or alternative model.

Path scores are comparable within a locus and are not probabilities.
