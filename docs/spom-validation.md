# S. pombe CNN-Mamba integration validation

Status on 2026-10-13: the three confirmed model/converter interface bugs are
fixed locally and covered by regression tests. After host GPU access was
enabled, real CNN-Mamba → UniAnn → AltAnn runs completed for all three nuclear
chromosomes on both strands. Saved native trace replay and recomputed Viterbi
outputs match byte for byte in every job. GFF3 structure and CDS phases pass,
but the biological audit identifies shared model limitations described below.
This establishes end-to-end execution, not universal biological validity.

## Inputs and provenance

Existing data were reused; no large dataset was downloaded or regenerated.
Paths below are relative to the project `DNA/` directory, containing `dev/`,
`data/` and `work/`. The assembly is RefSeq GCF_000002945.2 / ASM294v3, strain
972h-, assembled 2019-07-23, with PomBase annotation dated 2024-08-16.

| Input | Existing path | SHA-256 |
| --- | --- | --- |
| Genome | `data/basic_data/Spom/ref/genome.fna` | `fc3a90996a09d653880aea616e235c878a3224eb678f9d5f099b9d4821b29762` |
| Raw reference | `data/basic_data/Spom/ref/annotation.gff3` | `97ef0f5f89be5ac1bd117b0cf65abcadaf142b51f8c49be23515633157f5851c` |
| Frozen coding evaluation reference | `work/eviann_untrusted_cds_spom/reports/cds_loss_audit/protein_coding_reference.gff3` | `cab35be13c82e99fd266ae92ab89ef7efd09ec6d3ea0036851e85029b6f30031` |
| Available checkpoint | `dev/cnn-mamba/models/human_joint_frame_beta010_chr1held_epoch4_inference.pt` | `429904082463e9822ec83aa15cfd48b16cf0cdb6ff15ff4e7141316117f9c037` |
| Checkpoint calibration | `dev/cnn-mamba/results/human_joint_frame_beta010_v1/temperature_scaling/temperature.json` | `52633f955c6c0a5d16c3f78dce9dbd9c8414a51fcddd7a4ffe479ba6f44689cc` |

The genome/reference symlinks resolve into the existing
`GCF_000002945.2_dataset/ncbi_dataset/data/GCF_000002945.2/` archive.
Nuclear scope: NC_003424.3 (5,579,133 bases), NC_003423.3 (4,539,804), and
NC_003421.2 (2,452,883). Mitochondrial NC_088682.1 is excluded.

The available checkpoint is a **Human-trained** model. Its validation-fitted
temperature is 1.2288177735743935; this does not establish calibration or
prediction accuracy on S. pombe. cnn-mamba was inspected at
`809f74379efccf70ee7fd876dc95aabe5a4b0ea7`; local metadata changes are described
below. UniAnn was safely updated to `6fe3d46efacf1f900789a4cfa6d103c69dfe6ac0`,
then patched locally. Its native seven-state decoder source remains unchanged.

Existing S. pombe `AltAnn/example/` and work-directory scores use PSAURON
coding emissions and hibiki-kato/sitescore predictions. They are historical
compatibility data and were not substituted for Bowen's predictions.
No real four-state predictions were found in the initial data audit. They
were subsequently generated from the existing checkpoint and genome, with
explicit joint-head temperature 1.2288177735743935. No reference annotations
were supplied to model inference. An exact machine-readable manifest,
including resolved absolute paths and measured metrics, is saved under the ignored
`example/results/cnn-mamba-spom-validation/manifest.json`.

## Bugs and fixes

1. **Repeated calibration.** Bowen's scorer already supports multiclass
   `softmax(logits/T)` when `--calibration-json` or `--temperature` is supplied;
   its default is T=1 and calibration is not embedded in the checkpoint.
   UniAnn independently applied a binary T=7 transform to each probability,
   including already calibrated outputs. The local converter removes that
   transform and preserves the existing logarithmic scoring formulas, coding
   gate, stop penalties and two-decimal serialization. AltAnn performs no
   extra calibration. Whether any actual S. pombe prediction was calibrated
   was initially unknown; the real run now records the explicit temperature
   and calibration JSON in every joint-head scoring metadata file.

2. **Frame indexing.** Bowen class 0 is the first codon base; class 1 and 2
   are the next bases, advancing in transcript direction. GFF CDS phase `p`
   initializes the local label to `(3-p)%3`. For oriented zero-based position
   `n` and local class `i`, the UniAnn exon state is `(n-i)%3`. The old
   `(n-i-1)%3` maps it one frame behind. The fix follows independent biological
   checks against reference ATG/stop coordinates, phase continuity, and intron
   lengths on both strands, rather than changing labels by assumption.

3. **Header alignment.** Bowen exports an eight-column header followed by
   one row per base. UniAnn treated that header as a base. The corrected
   converter skips the exact exporter header, accepts headerless seven/eight
   columns, validates FASTA identifier, contiguous one-based positions,
   consistent strand, finite probabilities in [0,1], and exactly L rows.
   `predicted_state` is ignored.

4. **Stale emissions and hidden failures.** UniAnn reused an existing
   `out.ps.txt` even when the coding input changed; shell chains could hide
   failures behind the final success message. The wrapper now converts into a
   temporary file before replacing emissions, regenerates each run, and
   propagates converter/native failures. Tests verify that a failed conversion
   preserves the previous emission file and prevents native decoding.

5. **Published checkpoint metadata.** The site scorer indexed
   `checkpoint['optimizer_step']`, but the available inference checkpoint has
   epoch 4 and no optimizer step. That would fail while writing completion
   metadata after scoring. A small cnn-mamba helper records absent steps as
   null and preserves the epoch; three tests cover inference/training/missing
   progress cases. Neural inference code and logits are unchanged.

AltAnn's existing directory discovery was also fixed to prefer `out.ps.txt`
when a leftover PSAURON CSV is present, including with `-a`. Explicit legacy
raw inputs remain supported. No K-Viterbi, transition-rule or candidate-ranking
changes were necessary.

## Biological and downstream checks actually completed

The reference-coordinate audit checked 4,314 complete nuclear CDS models:
2,133 plus and 2,181 minus, comprising 8,672 CDS segments and 43,251 probes.
All pass ATG start, terminal stop, CDS length, internal-stop and GFF-phase
continuity checks. The old mapping disagrees at all 43,251 probes; `(n-i)%3`
agrees at all probes. The audit excludes 793 partial, 32 pseudo and 16
exception-bearing transcripts, plus 11 mitochondrial transcripts. These
biology-only exclusions do not redefine the frozen evaluation denominator.

Examples on chromosome I include plus-strand nhe1 (rna-NM_001018182.4), with
CDS 52117–52240 phase 0 and 52318–53600 phase 2, and minus-strand
rna-NM_001355965.2, whose transcript-ordered CDS intervals are 129172–129185
phase 0, 129103–129127 phase 1 and 128025–129059 phase 0. Both have ATG starts
and TAG terminal stops. The minus start is oriented position 5,449,948,
corresponding to E1, not the old E0. Seven real-DNA/GFF cases, covering all
initial frames on both strands, are frozen in UniAnn's regression fixture.
Probabilities in those tests are synthetic and annotation-derived.

The corrected synthetic-logit integration fixture uses the exact frozen
Bowen exporter TSV from the original contract investigation. The corrected
wrapper and freshly compiled UniAnn decoder ran on it; AltAnn consumed the
processed files. Tests check dense row alignment, frame evidence, unchanged
probability interpretation, saved-trace versus recomputed Viterbi equality,
literal float32 emission/transition score sums, alternative CDS start/stop
and phases, and reflected reverse output coordinates. The 12 corrected
fixture files were regenerated byte for byte. These checks establish score
compatibility, not real S. pombe candidate recall or final precision.

| Check | Before local fixes | After local fixes |
| --- | --- | --- |
| Rows from 1,700-base headered synthetic export | 1,701 | 1,700 |
| Reference frame probes agreeing with converter formula | 0 / 43,251 | 43,251 / 43,251 |
| Extra per-column T=7 transform | Applied | Removed |
| Changed coding input with existing emissions | Reused stale file | Regenerated |
| Published site-checkpoint metadata | Missing-key failure | Optional step, preserved epoch |

Passing tests: **80 AltAnn tests**, including the historical numerical
comparison with the original K-best executable and no skips; **12 UniAnn
converter/wrapper tests**; **3 cnn-mamba checkpoint-metadata tests**. Perl
syntax, shell syntax and diff whitespace checks pass. Previously completed
CPU calibration tests exercise Bowen's real temperature fitting utility.

## Commands used

From `dev/AltAnn`:

```sh
python3 scripts/audit_reference_frames.py \
  --genome ../../data/basic_data/Spom/ref/genome.fna \
  --reference ../../data/basic_data/Spom/ref/annotation.gff3 \
  --chromosome NC_003424.3 --chromosome NC_003423.3 --chromosome NC_003421.2 \
  --output example/results/cnn-mamba-spom-validation/reference_frame_audit.json \
  --fixture ../UniAnn/tests/data/spom_reference_frames.json

python3 tests/generate_corrected_fixture.py --uniann ../UniAnn \
  --gffread /home/hibiki/miniforge3/envs/eviann/bin/gffread \
  --output /tmp/bowen-corrected-regenerated

ALTANN_REFERENCE_BINARY=/home/hibiki/DNA/work/chimann/Dmel_k20_flank1000_chain/snapshot/bin/uniann \
  ctest --test-dir build --output-on-failure
```

From `dev/UniAnn`:

```sh
python3 -m unittest discover -s tests -p test_coding_frame_converter.py -v
perl -c scripts/convert_coding_frame.pl
bash -n scripts/uniann.sh
```

From `dev/cnn-mamba`:

```sh
PYTHONPATH=src python3 -m unittest discover -s tests -p test_checkpoint_provenance.py -v
/home/hibiki/miniforge3/envs/altann-spom/bin/python -c \
  'import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available(), torch.cuda.device_count())'
```

Initially the local Codex Bubblewrap sandbox overmounted `/dev` with a minimal
tmpfs, hiding NVIDIA devices. The kernel driver and host filesystem were
present; device lookup returned ENOENT rather than a permission-denied error.
After the user enabled full host access, `/dev/nvidia*` became visible,
`nvidia-smi` worked, and torch 2.10.0+cu128 detected the RTX 5070. Both a CUDA
matrix operation and the actual checkpoint's three output heads passed a GPU
smoke test. `pysam==0.23.3` was installed in the existing `altann-spom` environment;
no drivers or system security settings were modified by the agent.

Real run artifacts live in
`example/results/cnn-mamba-spom-validation/gpu-run-20261013/`. The run driver
`run.py` stages the corrected tracked scripts and a fresh native build, extracts
oriented nuclear FASTAs from the existing genome, invokes Bowen's existing
joint/site scorers and exporters, then the UniAnn wrapper and AltAnn. Every
command, working directory, exit code, log and duration is recorded in
`commands.jsonl`; per-job inputs, lengths and output counts are in
`run-manifest.json`. Parameters are K=10, flank=1,000, threads=4, GPU batch=4.
No legacy scores were used. Each model head was scored in 10 kb windows with
5 kb stride; joint-head calibration was explicitly enabled once upstream.

| Chromosome | Strand | Dense emission rows | Complete best-path transcripts | Alternative transcripts |
| --- | --- | ---: | ---: | ---: |
| NC_003424.3 | + | 5,579,133 | 764 | 1,776 |
| NC_003424.3 | - | 5,579,133 | 812 | 1,927 |
| NC_003423.3 | + | 4,539,804 | 644 | 1,497 |
| NC_003423.3 | - | 4,539,804 | 597 | 1,509 |
| NC_003421.2 | + | 2,452,883 | 316 | 673 |
| NC_003421.2 | - | 2,452,883 | 346 | 802 |

All six GFF/diagnostic TSV pairs are byte-identical between saved native trace
replay and recomputed Viterbi. The combined pool contains 11,663 transcripts,
including 3,479 complete best paths and 8,184 alternatives; five CDS chains
repeat across loci, so exact-CDS evaluation uses 11,658 unique chains.

UniAnn's wrapper retained two incomplete terminal models. Saved-trace analysis
identified only `sequence.fna.813` on chromosome I minus (terminal E2) and
`sequence.fna.598` on chromosome II minus (terminal I2) as incompatible with
AltAnn's documented complete-model contract. Separate `uniann.complete.gff3`
validation inputs exclude those IDs and their children; original wrapper
GFFs remain preserved. All other supplied complete models match the native
trace. This is explicit input preparation, not a new candidate filter.

## Evaluation conventions and measured results

The existing frozen coding reference contains 5,145 nuclear genes, 5,155
coding transcripts and 5,153 unique CDS chains. Follow the project's existing
`trmap -c '='` transcript comparison and `trmap -c '=' --strict-match` on
CDS-as-exons; gene recall counts a transcript OR strict complete-CDS match.
The exact CDS key is `(contig, strand, complete sorted CDS intervals)`.
No new coordinate tolerance is introduced. The full raw-reference evaluation
wrapper has noncoding/empty-ID denominator problems and is not a substitute
for the established coding-reference audit.

`scripts/evaluate_spom_predictions.py` reproduces the established convention,
including a regression against the historical audit's 3,801 matched genes,
4,127 unique query CDS and 3,675 exact reference CDS. No new tolerance or
reference-denominator filtering was introduced.

The controlled **before** comparison reuses the identical real calibrated
Bowen TSV and site scores with the original UniAnn `6fe3d46` converter,
including its frame offset and T=7 transform. Its header is explicitly removed
to allow a dense valid input; this is not a claim that the unmodified upstream
pipeline handles the header. `compare_original_converter.py` and
`original-converter-manifest.json` record this workaround, source hash and
parameters. Reference paths for this control are recomputed with AltAnn's
unchanged compatibility model; the corrected run additionally checks actual
native UniAnn traceback. There is no new model inference in the control.

| Metric | Original converter, header stripped | Corrected converter |
| --- | ---: | ---: |
| Best-path gene sensitivity | 15.45% (795 / 5,145) | 37.08% (1,908 / 5,145) |
| Best-path gene precision | 60.05% (795 / 1,324 loci) | 54.84% (1,908 / 3,479 loci) |
| Best-path CDS sensitivity | 15.12% | 35.88% |
| Best-path CDS precision | 58.84% | 53.15% |
| Best-path exact unique reference CDS | 779 | 1,849 |
| Candidate-pool gene recall | 18.56% (955 genes) | 45.56% (2,344 genes) |
| Candidate-pool exact CDS recall | 18.05% (930 chains) | 43.90% (2,262 chains) |
| Candidate-pool exact CDS precision | 25.17% | 19.40% |
| Alternative transcripts | 2,371 | 8,184 |

The corrected candidate pool retains all 3,479 best-path CDS chains and adds
8,179 unique chains, of which 413 exactly match additional reference CDS.
Gene precision of the candidate pool is 66.69% at source-locus level (2,320
matching source loci / 3,479); multiple candidates can recover multiple
reference genes from one source locus. This is not precision of a selected
final annotation. Final EviAnn selection was not run and its precision is
unmeasured. The before/after results jointly change frame indexing and scoring
calibration; they do not isolate either fix as an accuracy ablation. Recall
increases, while precision decreases with this Human checkpoint.

For reproduction from `dev/AltAnn`, the real run/control drivers and four
evaluation sets are preserved under the run directory. The evaluation command
uses `--reference` with the frozen reference above, repeated `--predictions`
for `baseline.gff3`, `candidates.gff3`, `original-converter-baseline.gff3` and
`original-converter-candidates.gff3`, `--candidate-label` for both pools,
all three `--chromosome` flags, and
`--trmap /home/hibiki/miniforge3/envs/eviann/bin/trmap`; it writes `evaluation.json`.
The biological audit command is:

```sh
python3 scripts/validate_prediction_cds.py \
  --genome ../../data/basic_data/Spom/ref/genome.fna \
  --gff example/results/cnn-mamba-spom-validation/gpu-run-20261013/candidates.gff3 \
  --output example/results/cnn-mamba-spom-validation/gpu-run-20261013/cds-audit.json
```

The audit deliberately returns 1 when biological findings exist while still
writing diagnostic JSON. `gffread -E` successfully loads all 11,663 transcripts.

## Remaining work and readiness

The real GFF has 61,095 feature rows and no structural, original-coordinate,
strand, phase, CDS-length, terminal-stop or serialized-score-delta failures.
11,546 of 11,663 transcripts pass all biological checks. The remaining 117
comprise 116 transcripts with internal stops assembled across splice junctions
(10 best paths, 106 alternatives) and one alternative with a non-ATG start.
There are 118 internal stop codons, all spanning junctions; none are wholly
within an exon. These models were reported, not silently removed or reranked.

The internal-stop limitation is shared with UniAnn: stop masking examines
contiguous genomic triplets, and its Viterbi states lack the partial-codon
history needed to identify stops assembled across introns. A native best-path
example joins chromosome I positions 1836783, 1836784 and 1838970 into TGA.
Real-DNA tests cover both strands in
`tests/data/spom_prediction_internal_stop.json`.

The non-ATG alternative is
`NC_003423.3.minus.kbest.locus235.k9.g1.t1`. Its first CDS exon is two bases:
the original oriented start motif ATG overlaps a donor GT at the G, so the
splice removes that G and the joined CDS starts ATC. UniAnn permits this
immediate first-exon intron transition. Site motifs and reverse mapping are
correct; the limitation is in the shared transition rules, not additional
AltAnn calibration or coordinate conversion. Fixing these shared model rules
would require changing the compatibility model beyond the three interface
fixes. A forced-path reproduction using the exact real DNA and observed site
scores produces the same two-base first exon in both native UniAnn and AltAnn;
its diagnostic emissions are synthetic and are not used in the evaluation.
`tests/data/spom_prediction_start_overlap.json` freezes the real-DNA case.
No such algorithm or ranking change was made.

For minus, each actual run scores an already reverse-complemented FASTA; the
exporter's strand flag alone would not transform sequence. Site tables label
each oriented job `+`, and AltAnn uses `--reverse` for original coordinates.

Select the model temperature explicitly and record scoring metadata. Use a
fresh scoring/output directory because Bowen's completion sentinel and TSV
export reuse completed outputs without checking a newly requested temperature.
The supplied Human calibration remains a transfer assumption, not verified
S. pombe calibration. Stage the corrected tracked UniAnn scripts with a fresh
native build: existing untracked `UniAnn/bin/` installation files were preserved
and may still contain the old converter.

ChimAnn's current
Nextflow AltAnn module still requests legacy raw PSAURON inputs and requires
a separate change to publish/pass these oriented standardized score files.
The corrected standalone pipeline has completed real end-to-end execution and
reference evaluation. It is usable for exploratory compatibility experiments.
Larger annotation-quality experiments need resolution of the shared biological
limitations, a suitable organism model/calibration, and the workflow migration.
The current run does not certify all exported CDSs as biologically valid.
