# Existing UniAnn file formats

AltAnn reads existing UniAnn inputs and optional outputs directly. No dedicated
export is required. Each job describes **one FASTA record in its scoring
orientation**, plus matching emission and motif scores. Minus-strand jobs must
already contain the reverse-complemented sequence and scores for that sequence.

## Sequence and retained scores

The FASTA identifier is the first token after `>`. DNA symbols are case
insensitive; IUPAC ambiguity symbols are accepted. Score positions refer to the
oriented sequence, before any conversion to chromosome coordinates.

`out.ps.txt` contains whitespace-separated rows:

```text
position  N  E0  E1  E2  I  [nucleotide]
0         0  -2  -2  -2  0  A
```

The first line above describes columns; it is not a literal file header.
Positions are zero based. Exactly one row is required for every FASTA position.
All five scores must be finite. The optional final column is one DNA symbol.
The intron score is shared by the internal `I0`, `I1`, and `I2` states.

Each of `out.gt.txt`, `out.ag.txt`, `out.atg.txt`, and `out.stop.txt` contains
`position score` pairs, also zero based. These are sparse tables of transformed
scores, not probabilities. Missing positions prohibit the corresponding
transition; empty files are allowed. Duplicate positions, extra columns, and
nonfinite scores are rejected. Retained tables allow blank lines and `#`
comments, but no un-commented header.

## Original PSAURON and site files

When retained emissions are absent, provide PSAURON's original
`psauron_score.csv` generated with `-a`. Its first four lines are the original
preamble, including the column header on the fourth line. Data rows require at
least 15 comma-separated columns. Column 1 matches the FASTA identifier;
columns 10 through 12 contain semicolon-separated forward-frame probabilities. Exactly
one row must match the FASTA record. Chromosome-sized fields are accepted.
AltAnn uses the bundled UniAnn Perl preprocessor to reproduce emissions.

When retained motif tables are absent, provide the original `sites*.tsv`.
Fields are whitespace separated:

| Column | Meaning |
| --- | --- |
| 1 | FASTA identifier; must match the job |
| 2 | One-based position |
| 3 | Orientation; only `+` rows are scored |
| 4 | `donor`, `acceptor`, `start`, or `stop` |
| 5 | Original auxiliary field; not used |
| 6 | Probability in the six-column form |
| 7 or later | Optional replacement probability; the last field is used |

An optional header beginning `chrom pos strand type` is accepted. Blank lines
and lines starting with `#` are ignored. Even a minus-strand job uses `+` site
rows relative to its already reverse-complemented FASTA.

Scores are `log(probability * multiplier + 1e-10) * factor`. A supplied UniAnn
log must contain one unambiguous `Multiplier is X, Factor is Y` pair. Otherwise,
the multiplier defaults to UniAnn's serialized `exp(1)`, and the factor is
estimated from the largest transformed donor score. For compatibility, factor
estimation uses column 7 when present, whereas final scoring uses the last
column. Use retained motif tables when an earlier run used custom preprocessing.
Log DP traces are not used to reconstruct numerical scores.

## Mapping and reference validation

A FASTA identifier matching `chromosome__start-end-total` describes a zero-based,
half-open segment in the original chromosome. Its length must equal `end-start`.
Other identifiers describe whole sequences unless `--seqid`, `--offset`, and
`--sequence-length` supply an explicit mapping. Reverse coordinates are
reflected within the segment, then shifted by its offset. GFF3 coordinates are
one based and inclusive. Sequence IDs and attributes are escaped for GFF3.

An optional `--best-gff` checks supplied transcript exon structures against the
recomputed reference before output filtering. A filtered subset, including an
empty result, is accepted. Nonempty models require transcript/mRNA and exon
features. The GFF may use local forward coordinates or mapped genomic
coordinates; all models must follow the same convention. It does not replace
the decoder's state path.

For batch input, a TSV manifest inventories existing files. Paths are relative
to the manifest directory; each row is a single oriented sequence. See the
[README](../README.md#both-strands-and-multiple-chromosomes) for columns and
examples. Output score attributes and sidecars are described in the
[output reference](../README.md#output).

## Runtime

Precompiled release packages contain the native decoder and its OpenMP runtime;
users do not need a compiler or a separate OpenMP installation. The frontend
requires Python 3.9 or later, with no additional Python packages. Perl is needed
only to regenerate emissions from PSAURON CSV files. Keep the archive's
`bin`, `share`, and runtime-library directories together when moving it.
