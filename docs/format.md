# UniAnn input formats

AltAnn reads processed UniAnn emission and transition scores. Each invocation
accepts **one FASTA record in one scoring orientation**. All supplied scores,
GFF records, and traceback entries must refer to that sequence in that
orientation. AltAnn does not reverse complement input files or generate scores
from PSAURON predictions.

## Sequence and processed scores

The FASTA identifier is the first token after `>`. DNA symbols are case
insensitive; IUPAC ambiguity symbols are accepted. Score positions refer to the
oriented input sequence, before any conversion to chromosome coordinates.

`out.ps.txt` contains whitespace-separated rows:

```text
position  N  E0  E1  E2  I  [nucleotide]
0         0  -2  -2  -2  0  A
```

The first line above describes columns; it is not a literal file header.
Positions are zero based. Exactly one row is required for every FASTA position.
All five scores must be finite. The optional final column is one DNA symbol.
The seven-state model is `N, E0, E1, E2, I0, I1, I2`: the intron score in the
file is shared by the three intron states.

Each of `out.gt.txt`, `out.ag.txt`, `out.atg.txt`, and `out.stop.txt` contains
`position score` pairs, also zero based. These sparse tables supply the
position-dependent transition scores used by UniAnn's seven-state model.
They contain transformed scores, not probabilities or a fixed 7-by-7 matrix.
Missing positions prohibit the corresponding transition; empty files are
allowed. Duplicate positions, extra columns, and nonfinite scores are rejected.
Tables allow blank lines and `#` comments, but no un-commented header.

All five score files are required. Use `--input DIRECTORY` to discover their
conventional names, or supply `--fasta`, `--emissions`, `--gt`, `--ag`, `--atg`,
and `--stop` explicitly. PSAURON CSV files and site probability tables are not
AltAnn inputs. Prepare the processed tables with UniAnn before decoding.

## Baseline path and scores

By default, supply both `--gff` (also accepted as `--best-gff`) and `--log`.
The log must contain the original UniAnn Viterbi `dp` and `bt` trace for the
same oriented sequence and score files. Each zero-based position has a `dp`
row followed by a `bt` row, each with exactly seven integer values. The initial
traceback row uses `-1`; subsequent predecessors are state indices `0` through
`6` in model order. A wrapper log containing only progress
messages or scaling parameters is insufficient. Traceback entries recover the
baseline path, including its noncoding regions and boundary history. GFF alone
cannot supply that path.

UniAnn truncates DP scores to integers. AltAnn uses the traceback and replays
its path with the supplied emission and transition tables to reconstruct
precise scores. It does not treat the printed integers as exact scores or
replace the supplied traceback with a newly optimized path. Incompatible or
incomplete traceback data is rejected. Rounded scores cannot independently
prove optimality: the supplied traceback is trusted, and its reconstructed
scores are checked against the log.

The GFF checks the baseline coding structure. It must describe the input
orientation, including for a reverse run. A filtered GFF is accepted when all
of its models match the recovered baseline. Incomplete terminal predictions
are not exported by AltAnn and must be excluded from a supplied reference GFF.

When a matching GFF or traceback log is unavailable, use `--rerun-viterbi`.
AltAnn then computes a new global baseline using its UniAnn-compatible model.
`--log` cannot be combined with this mode. A GFF is optional and, if supplied,
is checked against the computed baseline. The same local K-best search follows either
baseline source.

## Orientation and output coordinates

**`--reverse` requires a FASTA that you have already reverse complemented.**
Both the base order and the bases must be transformed (A with T, C with G).
Generate emissions and transitions for that transformed sequence; do not reuse
the forward scores. Any supplied UniAnn GFF and Viterbi log must come from the
same oriented input. AltAnn does not perform these input conversions.

Forward is the default. For a reverse run, prepare a reverse-complemented FASTA
and matching scores externally, then pass `--reverse`. Any supplied GFF and log
must also use that reverse input's coordinates. The flag changes output
coordinates and strand only. It does not reverse input arrays or infer their
orientation from filenames.

For an input sequence of length `L`, the one-based inclusive interval
`[start, end]` maps to `[L - end + 1, L - start + 1]` with strand `-`.
CDS phase is preserved. Forward output uses strand `+`.

A FASTA identifier matching `chromosome__start-end-total` describes a zero-based,
half-open segment in the original chromosome. Its length must equal `end-start`.
Other identifiers describe whole sequences unless `--seqid`, `--offset`, and
`--sequence-length` supply an explicit mapping. Reverse coordinates are
reflected within the segment, then shifted by its offset. GFF3 coordinates are
one based and inclusive. Sequence IDs and attributes are escaped for GFF3.

Each invocation accepts one sequence and one orientation. Decode additional
sequences or the opposite orientation in separate invocations. Output score
attributes and sidecars are described in the [output reference](../README.md#output).
The diagnostic TSV uses oriented input coordinates even with `--reverse`;
the main GFF uses mapped output coordinates.

## Model compatibility

AltAnn uses UniAnn's seven coarse states and retains per-path length and frame
history. That history does not introduce additional states. Length regulation
retains the current UniAnn constants: `MIN_INTRON = 40`, `MIN_EXON = 3`,
`MIN_INTER = 30`, and `MIN_SINGLE = 100`. These are fixed internally, as in
UniAnn; there are no command-line overrides. The original transition conditions
determine how each threshold applies.

Inputs from a modified UniAnn model may be incompatible with AltAnn's rules.
See [source provenance](../THIRD_PARTY.md) for the model lineage.

## Runtime

Precompiled release packages contain the native decoder and its OpenMP runtime;
users do not need a compiler or a separate OpenMP installation. The frontend
requires Python 3.9 or later, with no additional Python packages. Neither Perl
nor PSAURON is needed. Keep the archive's `bin`, `share`, and runtime-library
directories together when moving it.
