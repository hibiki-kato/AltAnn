# UniAnn input formats

AltAnn accepts one FASTA record with processed UniAnn emission and transition
scores, including scores derived from local coding-frame predictions. Legacy
PSAURON and site probabilities also remain supported. `-a` / `--all-prob`
accepts strand-tagged processed scores or legacy raw inputs.

## Sequence and processed scores (oriented inputs)

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
and `--stop` explicitly. This mode accepts existing UniAnn inputs and preserves
the saved GFF/traceback workflow.

Directory discovery prefers processed scores whenever `out.ps.txt` exists,
including with `-a` and a leftover `psauron_score.csv`. Explicit `--psauron`
or `--scores` selects legacy raw preprocessing. An incomplete processed score
set fails validation instead of silently substituting raw predictions.

## Current coding-frame conversion

UniAnn `uniann.sh -c coding.txt` uses `scripts/convert_coding_frame.pl` to
produce `out.ps.txt`. The locally corrected converter accepts Bowen's exact
eight-column header, or a headerless seven/eight-column table. Columns 4–6
are the three local coding-frame probabilities; column 7 is the noncoding
probability. It validates chromosome against the single FASTA identifier,
contiguous one-based positions, a consistent strand label, finite probabilities
in [0,1], and exactly one row per FASTA base. The optional eighth column is
ignored. Inputs must already use the matching sequence orientation.

For output position `n`, local probability column `i` (0, 1, 2) contributes
to absolute exon frame `(n - i) mod 3`. UniAnn performs this permutation
upstream: AltAnn reads `E0`, `E1`, and `E2` without permuting them again.
The exon states retain their absolute frame within the oriented sequence;
intron length determines the compatible frame on return to an exon.

The locally corrected converter performs no temperature scaling or probability
renormalization. With the supplied probabilities `q`, it uses
`log(exp(1) * q_noncode + 1e-10)` for N and normally I, and
`log(exp(1) * q_frame * (1 - q_noncode + 0.15) + 1e-10)` for E.
At the final base of a detected stop codon it sets the matching E score to
`-1000000` and I to `-2`. It writes two-decimal scores and the nucleotide.
AltAnn applies none of these transformations or temperature scaling again.
Model calibration belongs upstream; AltAnn consumes the resulting
standardized scores without needing model probabilities or logits.

The wrapper still transforms site probabilities into sparse zero-based
motif-start scores with `log(prob * mult + 1e-10) * factor`, where `factor`
is `round(1000 / max_positive_donor_log_score)`. Its raw site positions are
one based and only `+` rows are selected for an oriented job. `out.gt.txt`
and `out.ag.txt` identify the first GT/AG base; `out.atg.txt` and
`out.stop.txt` identify the first codon base. AltAnn consumes these scores
directly; `--mult` affects only its legacy raw-input preprocessing.

The current UniAnn wrapper disables its former `-a` helper. Score a reverse
complement separately upstream, then consume that job with `--reverse`.
Do not reuse plus-strand scores for the minus strand. The corrected wrapper
regenerates `out.ps.txt` atomically from the supplied coding table and propagates
conversion/native failures. Use separate job directories for each sequence
and orientation. Ordinary UniAnn tables have no strand column; the
combined format in the next section is AltAnn's additional input convention.

## cnn-mamba contract and discrepancies

The predictor is `bowencccc/cnn-mamba`, inspected at commit
`809f74379efccf70ee7fd876dc95aabe5a4b0ea7`. Its joint frame head emits four
logits in order `non_CDS, frame_0, frame_1, frame_2`. `P(CDS)` is
`1 - P(non_CDS)`; this is separate from the older three-class phase head,
which has no noncoding class. Use the four-state exporter with UniAnn.

`score_joint_frame_chromosome.py` applies `softmax(logits / T)` before writing
float32 arrays of shape `[sequence_length, 4]`. Calibration is selected with
`--calibration-json` or `--temperature`; omitting both uses `T=1`. The published
Human chr1-held-out model's fitted temperature is `1.2288177735743935`, from
validation logits. It applies only to the joint frame head, not the splice
or start/stop heads. This calibration is checkpoint-specific.
Use a fresh scoring directory when changing calibration: the scorer returns
early if its `.complete.json` exists, without validating the newly requested
temperature. The TSV exporter likewise reuses a completed output filename.

`export_joint_frame_tsv.py` writes a header and eight tab-separated columns:

```text
chrom position_1based strand frame0_probability frame1_probability frame2_probability noncoding_probability predicted_state
```

It reorders the four-class array into frame0/frame1/frame2/noncoding, without
another probability transform. Positions start at one. This agrees with the
probability column order expected by UniAnn. The unpatched `6fe3d46` converter
treats the header as a base, producing one extra emission row. The local fix
skips the exact exporter header; do not strip it for the corrected converter.
The final `predicted_state` column is ignored by UniAnn. Match chromosome,
orientation, row order, and sequence length before generating score files.

Bowen's `build_phase_track` defines local frames as codon offsets: a phase-0
CDS begins with labels `0,1,2,0,1,2`. On the minus strand those labels advance
in transcript direction. For an oriented base `n` with local codon offset `i`,
the corresponding absolute UniAnn frame is `(n - i) mod 3`. The unpatched
converter uses `(n - i - 1) mod 3`, shifting the coding evidence by one frame.
For example, a CDS beginning at zero-based position 60 has frame-0 evidence
mapped into E2 instead of the E0 state required by UniAnn's ATG transition.
This also affects the frame after a splice. The local fix uses `(n-i)%3`,
verified against complete S. pombe reference CDS on both strands, including
nontriplet introns and all initial frames. The earlier `coding-frame` fixture
uses synthetic probabilities arranged to match UniAnn's converter; it tests
score-file compatibility, not Bowen's actual label convention.

Unpatched UniAnn also applies its independent binary `T=7` transform after Bowen's
multiclass softmax, including when Bowen's fitted temperature has already
been applied. These two transformations are different, and the intermediate
four probabilities no longer sum to one after the binary transform. That
does not make the resulting emission scores unreadable, but it changes the
model's calibrated probabilities. The local fix removes that transform;
calibration is selected explicitly in Bowen's scorer. AltAnn neither recalibrates nor
relabels the saved scores, and its K-Viterbi model remains unchanged.

The chromosome scorer currently reads the FASTA in forward order and records
strand `+`. The TSV exporter's `--strand -` only changes a label; it does not
predict or reverse complement a sequence. For minus predictions, score an
already reverse-complemented FASTA upstream, use that job's processed scores,
and pass `--reverse` to AltAnn.

The real S. pombe run also confirms shared model limitations: contiguous
genomic stop masking does not catch stops assembled across splice junctions,
and first-exon intron entry can allow a donor to overlap the initial ATG.
Score compatibility and valid GFF3/CDS phases do not guarantee a biologically
valid CDS. See the validation report; AltAnn preserves the upstream transition
rules rather than adding candidate filtering or a different model.

The `bowen-contract` regression fixture runs the actual phase-label helper,
temperature utility, TSV exporter, and UniAnn converter on synthetic logits.
It freezes the original mismatches separately from tests that prove AltAnn's
literal score consumption and traceback replay. `bowen-corrected` replays the
same headered export through locally corrected UniAnn and checks native path
scores, alternatives, CDS phases, and reverse output coordinates. Neither
fixture contains real S. pombe model predictions. See the
[validation report](spom-validation.md) for the real reference-coordinate audit
and outstanding model inference/evaluation.

## Both-strand scores

Supply the original FASTA with `--fasta` (`-f`), a PSAURON CSV with `--psauron`
(`-p`), and site probabilities with `--scores` (`-s`). `-a` / `--all-prob`
prepares both orientations internally. With `--input DIRECTORY`, AltAnn
recognizes `psauron_score.csv` and `sites*.tsv`, `*_sites.tsv`, or `scores.txt`
alongside the FASTA when no processed emission file is present.
Explicit raw probability arguments cannot be mixed with processed score arguments.

Run PSAURON with its own `-a` option to retain all six frame probability arrays.
AltAnn reads the matching sequence row after the CSV's four-line preamble.
Columns 10-12 contain forward frame 0/1/2 probabilities; columns 13-15 contain
reverse-complement frame 0/1/2 probabilities. Each cell contains semicolon-separated
values. The reverse arrays are already in reverse-complement order and are used
without reversing their values or permuting frames. All probabilities must be
finite and between zero and one.

The site table contains six whitespace-separated columns:

```text
chrom  pos  strand  type      motif  prob
chr1   7    +       donor     GT     4.4e-05
chr1   900  -       donor     GT     3.0e-05
```

The header is optional; blank lines and `#` comments are accepted. `chrom` must
match the FASTA identifier. Rows for other sequences are ignored. `pos` is one
based in the original FASTA, within the supplied segment. It marks the first
motif base when read along the indicated strand, so a minus-strand motif uses
its rightmost genomic base. Types are `donor`, `acceptor`, `start`, and `stop`.
Rows may be interleaved and unsorted. An optional seventh probability column
replaces column six for scaling and scoring.

For a sequence of length `L`, minus site position `p` becomes `L - p + 1` in the
reverse complement. AltAnn applies UniAnn's emission preprocessing to the
oriented sequence and frame arrays, including stop-codon scores and two-decimal
rounding. It transforms site probabilities with `log(prob * mult + 1e-10)`.
`-m` / `--mult` ranges from 1 to 100 and defaults to exp(1). Each strand's factor
is inferred from its maximum positive transformed donor score. A strand without
a qualifying donor is skipped. Native score positions are zero based.

Both strands run through the existing seven-state native decoder independently.
Both global reference paths are recomputed, followed by local K-best decoding.
No traceback log or baseline GFF is required. `--log` and `--reverse` are rejected
with `-a`. An optional combined `--gff` must contain records in original genomic
coordinates; it is split by strand and checked against each computed baseline.
A skipped strand has an empty baseline and cannot match supplied GFF features.

The frontend merges predictions into one sorted GFF3 in original coordinates.
CDS phase, segment ownership, score attributes, and intron-chain selection are
unchanged. The diagnostic TSV has one header and a strand column; minus rows use
reverse-complement coordinates. JSON retains the original input paths for each
job. Without `-a`, raw inputs prepare only their supplied forward orientation.

Existing processed tables with a strand column also remain accepted with `-a`:
`position strand N E0 E1 E2 I [nucleotide]` for emissions and
`position strand score` for transitions. Their positions are zero based in the
original sequence; minus frame scores already refer to reverse-complement
frames. Numeric score text is copied unchanged into the oriented inputs.

## Baseline path and scores (oriented inputs)

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

## Orientation and output coordinates (oriented inputs)

**`--reverse` requires a FASTA that you have already reverse-complemented.**
Both the base order and the bases must be transformed (A to T, T to A, C to G, G to C).
For ambiguous IUPAC bases, exchange R with Y, K with M, B with V, and D with H;
S, W, and N remain unchanged. Apply the same mapping to lowercase bases.
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

Each invocation accepts one sequence. Both-strand mode decodes both
orientations; legacy mode decodes its prepared orientation. Output score
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

The native Viterbi model is unchanged between upstream `91477a6` and `6fe3d46`.
Inputs from a modified UniAnn model may be incompatible with AltAnn's rules.
See [source provenance](../THIRD_PARTY.md) for the model lineage.

## Runtime

Precompiled release packages contain the native decoder and its OpenMP runtime;
users do not need a compiler or a separate OpenMP installation. The frontend
requires Python 3.9 or later, with no additional Python packages. Neither Perl
nor PSAURON is needed. Keep the archive's `bin`, `share`, and runtime-library
directories together when moving it.
