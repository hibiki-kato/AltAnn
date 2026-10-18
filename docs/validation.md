# Validation

## Corrected UniAnn interface and S. pombe reference audit

Local UniAnn fixes accept Bowen's header, validate dense sequence alignment,
map local codon offsets with `(n-i)%3`, and remove temporary per-column `T=7`
scaling. The wrapper regenerates emissions and propagates failures. The
`bowen-corrected` fixture replays the original frozen synthetic-logit export
through these fixes and a fresh native decoder. All 80 AltAnn tests pass with
the historical numerical oracle enabled, with no skips. UniAnn's 12 converter
and wrapper tests and cnn-mamba's three checkpoint-metadata tests also pass.
The K-Viterbi implementation and candidate selection are unchanged.

An independent biological audit of the existing nuclear S. pombe reference
validates 4,314 complete CDS models (2,133 plus; 2,181 minus), 8,672 CDS
segments and 43,251 frame probes. Seven real-DNA/GFF cases are frozen in
`UniAnn/tests/data/spom_reference_frames.json`; test probabilities are derived
from annotations, so they are not model predictions.

After host GPU access was enabled and the missing pysam dependency installed,
real CNN-Mamba → UniAnn → AltAnn runs completed on all three nuclear chromosomes
in both orientations. Every saved-trace and recomputed-Viterbi GFF/TSV pair
matches byte for byte. The Human checkpoint gives 1,849 exact best-path CDS
matches and 2,262 in the candidate pool, with 8,184 alternative transcripts.
No legacy scores were substituted. GFF3 structure/phases pass, while 116 models
have splice-junction internal stops and one has an overlapping start/donor
motif producing a non-ATG CDS start; native reproductions demonstrate shared
UniAnn limitations. No candidate filtering or transition changes were made.
See [the dataset, metrics, commands, and remaining work](spom-validation.md).

## UniAnn coding-frame integration

UniAnn was updated to upstream `6fe3d46efacf1f900789a4cfa6d103c69dfe6ac0`.
Its native decoder source is unchanged from `91477a6`; the wrapper now uses
`convert_coding_frame.pl` for local coding-frame probabilities. The regression
fixture in `tests/data/coding-frame/` comes from that upstream converter,
site-score preprocessing, and native decoder, using synthetic probabilities.
It covers all three initial exon frames, intron lengths 60/61/62, fractional
emissions, stop penalties, alternative acceptors, and CDS phase 2 after splicing.
These are compatibility tests, not model-accuracy or calibration measurements.
The initial converter suite passed all 52 tests, including the historical
Dmel K-best comparison with `ALTANN_REFERENCE_BINARY` set to the original
patched binary.

A fresh installation staged in a temporary directory ran the updated UniAnn
wrapper end to end. Its five processed score files matched the regression
fixture byte for byte. AltAnn consumed the wrapper's saved `out.err` and
filtered annotation using directory discovery, with K=10 and flank=50,
producing three reference transcripts and four alternatives. ChimAnn's
`bin/merge_altann.py` accepted that output alongside the UniAnn annotation,
and gffread 0.12.8 accepted the merged GFF with seven transcripts and 13 CDS
features. The existing EviAnn integration consumes this merged GFF through `-c`;
EviAnn itself was not rerun for this small fixture.

Bowen's actual predictor is `cnn-mamba`, not the older `SpliceSiteCNN`
repository. It was checked out at `809f74379efccf70ee7fd876dc95aabe5a4b0ea7`.
Its CNN/bidirectional-Mamba backbone feeds splice and start/stop heads, plus
a dilation-3 frame branch with a four-class joint frame head. The published
checkpoint was inspected on CPU: it has eight blocks, feature width 192,
kernel size 7, 10 kb windows, 5 kb stride, and a `[4,192]` joint-head weight
matrix. It contains no embedded temperature; calibration is a separate JSON.
It implements validation-fitted multiclass logit temperature scaling. The published Human
model's temperature is `1.2288177735743935`; inference uses `T=1` unless a
temperature or calibration JSON is supplied explicitly.

Comparing its actual label and export implementation with UniAnn identified
three upstream discrepancies: the exported TSV header is treated as a base;
the converter's `(n-i-1)%3` permutation is one frame off relative to Bowen's
first-codon-base-is-0 labels; and the converter's independent `T=7` probability
transform runs even after model calibration. The original synthetic converter
fixture deliberately follows the converter's permutation and cannot establish
that Bowen's labels are interpreted correctly.

The additional `tests/data/bowen-contract/` fixture uses synthetic logits with
the actual phase-label function, NumPy temperature utility, and TSV exporter,
followed by the unmodified UniAnn converter and native decoder. Its diagnostic
tests record the original header, frame, and temperature discrepancies; its
decoding regressions check that AltAnn preserves those saved score files and
reproduces the native traceback exactly. Runtime tests remain independent of
cnn-mamba, NumPy, pandas, and PyTorch. Regeneration uses those packages on CPU.
CUDA was unavailable during this initial fixture investigation; the later
real GPU run is reported above.
AltAnn still requires no direct model dependency and adds no calibration.
All 59 AltAnn tests passed, including the historical reference comparison,
with no skips. Regeneration reproduced all 14 Bowen contract artifacts byte
for byte. Two existing cnn-mamba calibration tests (ambiguous-label masking
and temperature fitting/argmax preservation) also passed on CPU, executing
their unchanged methods separately from the GPU-dependent model suite.

The current UniAnn wrapper disables its former both-strand helper. Its existing
installed `bin/` files were preserved during the repository update; the test
installation used the updated tracked scripts and a fresh native build.
ChimAnn's current `modules/altann.nf` still invokes legacy raw PSAURON mode.
Migrating that workflow to publish and pass the new oriented score files is
separate from this standalone AltAnn compatibility change.

## Version 0.3.0 both-strand validation

The chromosome III tutorial reads a six-frame PSAURON CSV and a combined site
probability table with the original FASTA. With K=10, flank=1,000, and four
threads, its GFF3 and TSV are byte-identical to the two independently oriented
processed-score runs after merging. The output contains 13,688 feature rows and
2,222 transcripts: 787 references and 1,435 alternatives, with 1,114 plus and
1,108 minus transcripts. The example probabilities were combined from the
previously scored orientations without rerunning their scoring models.

A separate 60,000-base sequence was scored directly with PSAURON 1.1.3 `-a`.
All five processed score files on each strand match UniAnn's original Perl and
shell preprocessing byte for byte. The both-strand raw AltAnn run matches two
independent processed-score runs in all GFF fields and TSV scores. Randomized
emission inputs, transition probabilities, and donor-factor rounding boundaries
also match the original preprocessing numerically.

Regression tests cover raw six/seven-column site probabilities and existing
processed inputs, all sequence-length remainders, unsorted rows, distinct reverse
frame probabilities, donorless skips, nonzero reverse CDS phase, combined GFF
validation, segment mapping, invalid-input preservation, and deterministic
threaded output. Installed and extracted release packages exercise raw `-a`
inputs alongside the existing processed-score mode. The C++ decoder and its
seven-state model are unchanged.

## Version 0.2.0 saved traceback validation

The current tutorial inputs for S. pombe chromosome III were decompressed into
a fresh directory. The original seven-state UniAnn executable at upstream
commit `91477a69e1a949fed082c4662b348d9a7a91ca2b` produced a GFF and full DP/BT
log separately for the forward and reverse-complemented inputs.

For each direction, AltAnn decoded those original logs with four threads and
also ran `--rerun-viterbi` with one thread, using K=10 and flank=1,000. The final
GFF and diagnostic TSV files were byte-identical between the two modes. The
reverse run used `--reverse`; input scores and traceback remained unchanged.

| Direction | Complete references | Alternatives | Total transcripts |
| --- | ---: | ---: | ---: |
| Forward | 385 | 729 | 1,114 |
| Reverse | 402 | 706 | 1,108 |

The original logs are about 220 MB per direction. These logs, original GFFs,
and comparison results remain ignored under `example/work/` and
`example/results/`; they are not part of the distributed example.

The local suite passes all 25 tests, including the optional original snapshot
comparison. The saved-log regression uses checked-in files generated by the
original UniAnn binary, checks both output orientations, and rejects corrupted
scores, predecessor changes, missing rows, and incompatible CLI arguments.

The following sections retain version 0.1.1 measurements. Their historical
whole-genome commands are not the current chromosome III tutorial.

## Numerical compatibility

The Dmel reference is the saved `Dmel_k20_flank1000_chain` coarse-state decoder
and its existing chromosome-wide outputs. The metadata-expanded implementation
is not used as an oracle.

Synthetic integration checks cover independent expected exon coordinates,
multiple splice choices and intron retention, K limits, score deltas, byte-identical
serial/parallel output, both strands and nonzero CDS phase, segment ownership,
manifest ordering, and malformed inputs. The optional snapshot test compares
raw GFF and every TSV field against the original executable.

On Dmel chromosome 4 (1,348,131 bases), K=20 and flank=1,000:

- Raw positive-strand GFF, reference GFF, and TSV were byte-identical between
  AltAnn with 1 and 4 threads, a fresh original decoder run, and archived results.
- Complete frontend runs regenerated scores from the original FASTA, PSAURON,
  site tables, and logs for both strands. With the short-CDS filter enabled,
  final GFF feature records matched the archived output exactly after sorting:
  7,183 records on plus and 7,294 on minus (357 and 478 transcripts).
- Genome-wide equality across the other six chromosomes has not been measured.

## S. pombe example

All three nuclear chromosomes were scored in both orientations for the
original version 0.1.1 example. The original seven-state UniAnn
binary at upstream commit `91477a69e1a949fed082c4662b348d9a7a91ca2b` was run
on all six prepared inputs. AltAnn 0.1.1, K=10, flank=1,000, and four threads
reproduced every complete UniAnn CDS chain after reverse-coordinate mapping:

| Chromosome | Plus references | Minus references |
| --- | ---: | ---: |
| NC_003424.3 | 974 | 1,060 |
| NC_003423.3 | 840 | 807 |
| NC_003421.2 | 385 | 402 |

UniAnn also emitted one terminal partial model on chromosome I minus.
AltAnn excludes this incomplete model as documented. The final output contains
4,468 complete references and 8,226 alternative transcripts, after intron-chain
deduplication. The full frontend run from the compressed example inputs took
98.3 seconds on this development machine in one measurement. This includes
score preprocessing and GFF conversion, but excludes archive decompression.

This example exposed a flank extending into an omitted terminal partial model.
The boundary regression now checks that the last complete reference is kept,
the partial coding segment is excluded, and threads 1 and 4 produce identical
output. Separate checks cover original UniAnn CDS-only GFF validation in both
coordinate conventions. The version 0.1.1 local suite passed all 21 tests when the original
snapshot executable was supplied.

## Dmel timing

Single-run native-decoder measurements on this development Linux amd64 machine,
using the same prepared chromosome-4 plus inputs and `/usr/bin/time`:

| Decoder | Threads | Elapsed seconds | Peak RSS (KiB) |
| --- | ---: | ---: | ---: |
| Original snapshot | 1 | 11.93 | 532,728 |
| AltAnn | 1 | 4.13 | 158,548 |
| AltAnn | 4 | 2.08 | 214,972 |

These runs have 26 reference loci. They exclude preprocessing and Python GFF
conversion and are not a whole-genome speed guarantee. Local traceback storage
grows with K and locus length, and simultaneous loci increase peak memory.

## Platform validation

The CI matrix builds, tests, installs, and packages Linux amd64/arm64 and macOS
amd64/arm64 natively. Runtime bundling is checked separately from source builds.
See the repository's Actions results for the current commit's platform status.
Runner labels follow the [GitHub-hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
