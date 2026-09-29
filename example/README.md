# S. pombe example

This example runs AltAnn on all three nuclear chromosomes of
*Schizosaccharomyces pombe*, on both strands. The compressed inputs include
the genome, PSAURON coding scores, and predicted splice, start, and stop site
probabilities. The mitochondrial chromosome is excluded from decoding.

The compressed inputs total 43.3 MB. File checksums are listed in
`SHA256SUMS`, and software versions and model hashes are recorded in
`provenance.json`. Extracted inputs occupy about 335 MB, with additional
temporary space used while decoding.

Use AltAnn 0.1.1 or later, Python 3.9 or later, and Perl with either a built checkout or an
extracted release package. Perl converts the supplied PSAURON CSV files.
See the [installation instructions](../README.md#install-a-release) for the
release packages and [build instructions](../README.md#build-from-source-developers)
for a source checkout. The supplied scores are ready to use; no GPU, model
training, or scoring software installation is needed.

Run these commands from the repository root after building AltAnn:

```sh
python3 example/prepare.py
bin/altann decode --manifest example/work/inputs.tsv \
  --k 10 --threads 4 --output example/results/spom.gff3
```

With a release package, use its `bin/altann` path in the second command.
The preparation step decompresses the score tables and writes a FASTA for
each chromosome and strand. Minus-strand FASTAs are reverse complemented;
their supplied scores use the same orientation. `inputs.tsv` lists all six
jobs. AltAnn converts the final GFF coordinates back to the original genome.

`--k 10` retains up to ten paths per state during local decoding. Use
`--k 5` for fewer candidates. K does not specify the number of exported
transcripts: the UniAnn best is retained as rank 0, duplicate intron chains
are removed, and a path can contain more than one gene. `--threads 4` sets
the number of workers for local searches; `--threads auto` uses the available
CPU count.

The run writes three files:

| File | Contents |
| --- | --- |
| `results/spom.gff3` | Transcripts, exons, and CDS features in genome coordinates |
| `results/spom.gff3.tsv` | Candidate diagnostics in oriented input coordinates, including candidates removed during export |
| `results/spom.gff3.json` | Input paths, parameters, scaling, counts, and elapsed times |

These paths are relative to `example/`. Both `work/` and `results/` are
ignored by Git.

With the commands above, AltAnn 0.1.1 exports 12,694 transcripts: 4,468 complete
UniAnn references and 8,226 alternatives. One incomplete chromosome-end model
is excluded. Changing K or filtering options changes the candidate count.

GFF column 6 is `.`. Transcript attributes contain `kbest_score`,
`kbest_reference_score`, `kbest_delta`, `kbest_rank`, and `kbest_origin`.
The delta is the candidate path score minus the UniAnn best score over the
same local interval. Scores describe whole paths and should be compared
within the same locus; they are not probabilities. See the
[output reference](../README.md#output) for all attributes.

The genome is the NCBI RefSeq assembly `GCF_000002945.2` (ASM294v3).
The decoded chromosomes are `NC_003424.3`, `NC_003423.3`, and `NC_003421.2`.
Coding probabilities were generated with PSAURON 1.1.3. Site probabilities
were generated with the saved ConvMamba model trained for S. pombe in the
ChimAnn example. The genome is NCBI reference data, and the score tables are
derived model predictions. AltAnn's GPLv3 license and UniAnn attribution
cover the software; see [source credits](../THIRD_PARTY.md) for its origins.
