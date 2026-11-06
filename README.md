# AltAnn

Alternative gene annotations through local K-best Viterbi decoding.

AltAnn reads standardized UniAnn scores and produces scored alternative gene models
in GFF3. A C++17 core handles dynamic programming and parallel local searches;
a Python standard-library frontend prepares probabilities or reads processed scores and exports
annotations in original genome coordinates. Each invocation decodes one FASTA
record, with strand-tagged processed inputs supported by `-a` / `--all-prob`. No UniAnn
executable is required.

AltAnn builds on **[UniAnn by Aleksey Zimin](https://github.com/alekseyzimin/UniAnn)**
and inherits its **GPLv3** license. Original author credits are retained;
see [AUTHORS.md](AUTHORS.md) and [source provenance](THIRD_PARTY.md).

Try the [S. pombe example](example/README.md) to decode chromosome III using
supplied score files in one both-strand run.

## Install a release

Release packages target Linux/macOS on amd64/arm64. Download the matching
archive from [Releases](https://github.com/hibiki-kato/AltAnn/releases), extract
it, and run its `bin/altann` launcher. The native decoder and its OpenMP runtime
are included: **users do not need a compiler or a separate OpenMP installation**.

```sh
tar -xzf altann-linux-amd64.tar.gz
altann-linux-amd64/bin/altann decode --input /path/to/uniann-job \
  --rerun-viterbi --k 10 --threads auto --output alternatives.gff3
```

Python 3.9 or later is needed for the frontend, with no additional Python
packages. Perl and PSAURON are not required. Release archives are directory
packages containing a native executable and a Python frontend. Linux release builds use Ubuntu 22.04 (glibc 2.35);
macOS packages are built on macOS 15. Older systems have not been validated.

## Build from source (developers)

Targets: **Linux and macOS, each on amd64 (x86_64) and arm64 (aarch64)**.
Requirements: a C++17 compiler with OpenMP, CMake, and Python 3.9 or later.

```sh
git clone https://github.com/hibiki-kato/AltAnn.git
cd AltAnn
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
bin/altann --help
```

On macOS, install the build tools and OpenMP runtime first:

```sh
xcode-select --install
brew install cmake libomp python
```

Both Apple Silicon and Intel Homebrew prefixes are detected. For a custom
OpenMP installation, configure with `-DCMAKE_PREFIX_PATH=/path/to/libomp`.
Linux users need their distribution's C++ compiler and OpenMP development
package. The optional Makefile is a shortcut for GCC on Linux (`make -j4`).

Optional installation into a user-owned prefix:

```sh
cmake --install build --prefix "$HOME/.local"
```

The source build avoids `-march=native` so it does not deliberately restrict
the executable to the build machine's CPU. CI builds and tests all four target
OS/architecture combinations. Main-branch builds produce downloadable package
artifacts. After the tested version is merged into main, CI automatically tags
and publishes a new version as a GitHub Release. Existing versions are not
overwritten; see [maintenance](docs/maintenance.md).

## Current UniAnn coding-frame inputs

UniAnn upstream `6fe3d46` replaces PSAURON preprocessing with
`scripts/convert_coding_frame.pl`, using Bowen's local coding-frame probabilities.
Its output is still `out.ps.txt`; the five emission columns and seven-state
Viterbi model are unchanged. AltAnn consumes these processed emissions and
the four `out.*` transition tables directly, without additional calibration or
score transformation. Neither cnn-mamba nor a global-frame prediction
file is required by AltAnn.

The actual predictor is [Bowen's cnn-mamba](https://github.com/bowencccc/cnn-mamba).
The unpatched UniAnn `6fe3d46` interface has three discrepancies: header handling,
a one-frame indexing offset, and repeated temperature scaling. Local UniAnn
fixes accept Bowen's header, map local class `i` at position `n` to `(n-i)%3`,
and remove the temporary `T=7` transform. Real S. pombe reference CDS coordinates
confirm the frame convention on both strands. AltAnn consumes the corrected
scores without changing its K-Viterbi algorithm. Real inference and reference
evaluation completed on all three S. pombe nuclear chromosomes in both
orientations after host GPU access was enabled. The audit also identifies
shared UniAnn limitations involving splice-junction stops and overlapping
start/donor motifs; see the [S. pombe validation report](docs/spom-validation.md).
See the [verified model contract](docs/format.md#cnn-mamba-contract-and-discrepancies).

```sh
bin/altann decode --input /path/to/uniann-job --rerun-viterbi \
  --k 10 --threads auto --output alternatives.gff3
```

Use the FASTA from that exact UniAnn job. To reuse its saved reference path,
supply `--gff genome.fa.uniann.gff --log out.err` instead of
`--rerun-viterbi`. For a separately scored reverse-complemented FASTA, add
`--reverse` to map the output back to the original chromosome. Run each
chromosome and orientation separately with the current UniAnn outputs.
UniAnn's latest wrapper disables its former both-strand helper; its ordinary
score tables have no strand column. AltAnn's `-a` requires combined
strand-tagged processed tables or legacy PSAURON inputs as described below.

Directory discovery prefers `out.ps.txt` even if a legacy `psauron_score.csv`
is also present. Raw PSAURON preprocessing remains available through explicit
`--psauron` and `--scores` arguments. See the
[coding-frame score contract](docs/format.md#current-coding-frame-conversion)
for upstream calibration and coordinate details.

## Legacy PSAURON inputs on both strands

Supply the original FASTA, a six-frame PSAURON CSV generated with PSAURON
`-a`, and a site probability table containing both `+` and `-` rows:

```sh
bin/altann decode -f chromosome.fna -p psauron_score.csv -s sites.tsv -a \
  --k 10 --threads auto --output alternatives.gff3
```

`-a` (`--all-prob`) reverse complements the FASTA and prepares emission and
transition scores internally using UniAnn's preprocessing rules. It computes
an independent global reference path for each, and writes one sorted GFF3.
The diagnostic TSV contains both strands; JSON records both jobs and the
original input paths. No baseline log or GFF is required. An optional combined
`--gff` is checked against the computed references; `--log` and `--reverse`
cannot be combined with this mode.

The site table uses `chrom pos strand type motif prob` columns, with an optional
seventh rescaled probability column. Positions are one based in the original
FASTA and mark the first motif base when read along the indicated strand.
Sequence identifiers must match the FASTA. `-m` (`--mult`) controls probability
rescaling as in UniAnn; the default is exp(1). Strands without a donor with a
positive rescaled score are skipped. See the
[detailed input contract](docs/format.md#both-strand-scores).

## Decode an existing UniAnn run

```sh
bin/altann decode --input /path/to/uniann-job \
  --gff best.gff --log uniann.log \
  --k 10 --flank 1000 --threads auto --output alternatives.gff3
```

The input contains exactly one FASTA record and its processed UniAnn scores.
All files must describe the same sequence in the same orientation.

| Input | Files |
| --- | --- |
| Sequence | One `.fa`, `.fasta`, or `.fna` file |
| Emissions | `out.ps.txt` |
| Transition scores | `out.gt.txt`, `out.ag.txt`, `out.atg.txt`, `out.stop.txt` |
| Baseline annotation | UniAnn GFF, supplied with `--gff` (`--best-gff` is an alias) |
| Baseline traceback | UniAnn Viterbi log, supplied with `--log` |

Ambiguous file selection fails with an explanation. Explicit paths override
score discovery: `--fasta`, `--emissions`, `--gt`, `--ag`, `--atg`, and `--stop`.
This mode consumes UniAnn's already processed emission and transition scores.
Use the both-strand command above to supply PSAURON and site probabilities directly.

By default, GFF and log are required. AltAnn follows the log's traceback to
recover the baseline path, then replays that path using the supplied scores.
UniAnn's printed DP scores are truncated to integers, so they cannot provide
exact path scores. The GFF is checked against the recovered baseline. Filtered
annotations are accepted as matching subsets; incomplete terminal models must
be excluded from a supplied reference GFF;
see [input formats](docs/format.md).

If the GFF or traceback log is unavailable, explicitly request a new global
Viterbi pass:

```sh
bin/altann decode --fasta chromosome.fna --emissions out.ps.txt \
  --gt out.gt.txt --ag out.ag.txt --atg out.atg.txt --stop out.stop.txt \
  --rerun-viterbi --threads 8 --output alternatives.gff3
```

In this mode, `--log` is not accepted and the GFF is optional. A supplied GFF still
checks the computed baseline. Both modes use the same seven-state model and
local K-best search.

## Existing oriented inputs

**Before using `--reverse`, you must supply a reverse-complemented FASTA.**
Reverse complement means reversing the sequence and replacing each base with
its complement (A to T, T to A, C to G, G to C). Reversing the base order alone is insufficient.
This flag only maps the output; the both-strand mode above prepares both
orientations from the original sequence.

In this existing mode, one invocation processes one orientation. To annotate
the opposite strand,
first reverse complement the FASTA outside AltAnn and produce UniAnn scores
for that sequence. Supply those files in a separate invocation with `--reverse`:

```sh
bin/altann decode --input /path/to/reverse-uniann-job \
  --gff reverse.best.gff --log reverse.uniann.log --reverse \
  --threads 8 --output reverse.alternatives.gff3
```

The reverse FASTA, emissions, transition scores, GFF, and log must all use the
reverse input's coordinates. AltAnn trusts that these inputs have been prepared
consistently. `--reverse` changes the output mapping only: annotations are
reported on the original sequence with strand `-`. It does not transform
input sequence, scores, or traceback. Without the flag, output uses strand `+`.
Orientation is never inferred from filenames.

For a complete sequence of length `L`, an input interval `[start, end]` becomes
`[L - end + 1, L - start + 1]` in the output GFF. These coordinates are one based
and inclusive. CDS phase is preserved. Use `--seqid` if the reverse FASTA has a
different identifier from the original chromosome.

Headers of the form `chromosome__start-end-total` carry zero-based half-open
segment coordinates. Other headers use the complete sequence by default;
`--seqid`, `--offset`, and `--sequence-length` specify an explicit mapping.
Reverse coordinates are reflected within the supplied segment, then shifted
by its offset. For segmented input, midpoint ownership and internal-edge
exclusion use `--segment-overlap 4000000` and `--segment-margin 20000`. Set these
to the segmentation settings used when preparing the inputs.

Run additional sequences separately. The both-strand mode combines both
orientations of the same sequence in one output.

## Model and selection

The compatibility model has seven coarse states:
`N, E0, E1, E2, I0, I1, I2`. The emission table has five scores per position:
`N, E0, E1, E2, I`; the intron score is shared by all three intron frames.
An optional final nucleotide column is accepted. Positions in score tables
are zero-based; sparse transition files contain `position score` pairs. These
position-dependent scores work with UniAnn's transition rules; they are not a
fixed 7-by-7 transition matrix.

Length regulation uses the current UniAnn constants: `MIN_INTRON = 40`,
`MIN_EXON = 3`, `MIN_INTER = 30`, and `MIN_SINGLE = 100`. UniAnn fixes these
values internally, so AltAnn does too. They are model constants, not additional
states or command-line parameters. The original transition predicates determine
how the thresholds apply.

AltAnn retains K paths per coarse state (`--k`, default **10**; use `--k 5`
for a smaller candidate set). Length and frame history accompany
each path to enforce UniAnn's transition rules, but **do not expand the state
space**. This preserves the Dmel compatibility search; it is not a guarantee
of the mathematical global top K over all history-dependent paths.

The default search includes each reference gene plus at most 1,000 intergenic
bases on each side, clipped by neighboring coding states and sequence ends.
This also excludes incomplete terminal genes from a neighboring locus's flank. Both local
end states are N. Boundary history comes from the global N prefix. A start
codon may not extend before the local window. Complete genes are exported;
partial terminal predictions are not completed by this local search.

An alternative is retained only if its CDS shares at least one genomic base
with the CDS of its own global-best reference transcript, on the same sequence
and strand. Overlapping start-to-stop spans alone are insufficient: independent
ORFs entirely inside a reference intron or in the intergenic flanks are excluded.
This overlap rule always retains the global-best reference. A local path can still
contain other genes; only overlapping alternatives are ranked and exported,
while their scores remain scores of the complete local path. Diagnostic
`genes_in_path` and `gene_in_path` describe that original path.

The diagnostic TSV column `splice_site_jaccard` compares each candidate with
its global-best reference: donor and acceptor positions are compared as
separate sets, |shared| / |union|. Start and stop codons are ignored; two
intronless models score 1, and an intronless model against a spliced one
scores 0. It is reported only; no candidate is filtered on it. The final
GFF3 carries it as `kbest_splice_jaccard`, rounded to six significant
digits (1 for references).

Reference transcripts are retained with rank 0. Candidates use the original
local decoder ranks, up to K. Intron chains are deduplicated within each locus;
the reference takes precedence over alternatives sharing its chain. Empty
chains (single-exon models) also deduplicate within the locus. K is the number
of retained local paths, not a promise of K exported isoforms. A path may
contain several genes; disjoint genes are excluded before annotation ranking.
Selection happens before optional output filtering.

`--raw` turns off that selection. Every gene of every decoded path is exported,
including genes that share no CDS base with the reference, repeated annotations
and paths matching the reference, and `rank` is the decoder's path index. Paths
without a complete gene write nothing but still consume a rank. Genes without
CDS overlap have `classification` `no_cds_overlap` and `splice_site_jaccard` 0.
Segment ownership still applies, so overlapping segments do not repeat a locus.
`--raw` cannot be combined with `--short-cds-filter`.

`--short-cds-filter` reproduces the Dmel export rule: omit transcripts with at
most two exons and at most 200 coding bases. It is disabled by default so this
filter is an explicit choice.

The model follows the coarse-state UniAnn snapshot documented in
[THIRD_PARTY.md](THIRD_PARTY.md). Its native state definitions and transition
rules remain unchanged in upstream `6fe3d46`; the new coding-frame converter
feeds that same model. Customized scoring rules must be checked before use.

## Output

`alternatives.gff3` contains transcript/exon/CDS records. GFF column 6 is `.`;
scores are attributes, matching the Dmel whole-genome export convention:

| Attribute | Meaning |
| --- | --- |
| `kbest_score` | Score of the complete local path |
| `kbest_reference_score` | Baseline path replayed over the same local interval |
| `kbest_delta` | Stored candidate score minus stored reference score |
| `kbest_rank` | 0 for reference, 1..K for candidates |
| `kbest_origin` | `uniann_best` or `local_kbest` |
| `kbest_reference` | `uniann_best` |
| `kbest_locus` | Unique source segment/strand/locus identifier |
| `kbest_gene_in_path` | Gene identifier within a ranked path |
| `kbest_splice_jaccard` | Donor/acceptor Jaccard against the global-best reference |

Scores use 17 significant digits. Deltas are computed with decimal arithmetic
from the serialized scores. These are path scores, not probabilities or
independently scored isoforms. Compare them within a locus. Several genes in
one path share its score. Positive deltas against the legacy global reference
are possible.

Two sidecars accompany the GFF:

- `.gff3.tsv`: native candidate diagnostics, with segment and strand columns.
  Coordinates in this report are in the oriented input sequence. It includes
  candidates removed by export selection; rank-0 references are in the GFF.
- `.gff3.json`: input paths, parameters, counts, and elapsed times.

The main GFF is atomically replaced after successful decoding and conversion.
Input files are never modified.

## Score supplied models

`altann score` assigns a model score to existing CDS models, for example
EviAnn isoforms, instead of searching for alternatives:

```sh
bin/altann score --input /path/to/uniann-job --rerun-viterbi \
  --models eviann.gff3 --threads auto --output eviann.scores.tsv
```

It accepts the same processed inputs, `--log`/`--rerun-viterbi`, `--reverse`
and segment options as `decode`. Transcript CDS chains are read from GFF3
(`Parent=`) or GTF (`transcript_id`); only models on this sequence and
orientation are scored, and the JSON sidecar counts the others. Each model
replaces the global-best path over its span, extended to every global coding
run it touches, and both paths are replayed over the decoder's N-bounded
window. `delta` is therefore the model score minus the global-best score on the
same interval: 0 for a global-best gene, and identical to `kbest_delta` for a
single-gene K-best candidate. Models the seven-state rules cannot represent are
reported as `forbidden`, with the first failing transition (`start`, `donor`,
`acceptor`, `stop` or `internal_stop`) and its genomic position; malformed
models are `unscorable`. Non-canonical donors and stops spanning a splice
junction are common reasons, so `forbidden` means "outside the model", not
"false".

## Performance and development

`--threads auto` uses CPUs available to the process. `--threads N` provides
explicit control. Memory grows slowly with sequence length, about 10 MB per Mb,
and with the number and length of loci decoded concurrently; use a smaller
count on memory-constrained machines. Auto currently selects CPUs, not a memory
budget. A 235 Mb test sequence peaked at 2.4 GB with K = 10, so chromosomes need
no splitting for AltAnn itself.

A locus's traceback takes `28 x K` bytes per window base. `--traceback-memory`
(default 4096 MB) is shared by all loci decoding at the same time, however many
threads run. A locus that does not fit moves its traceback through a temporary
file beside the output, using background writes and read-ahead, and holds at most
32 MB in memory meanwhile. The file is unlinked when created, so it never
outlives the process, but it needs free disk space for every locus spilling at
once. This choice does not affect the output.

The global pass walks the emission file forward with rolling scores and a compact
traceback. The local pass reads only each locus's window from that file,
parallelizes independent loci with OpenMP and writes results in fixed order.
The Python frontend sorts the exported annotations. A sequence's candidate
annotations are held in memory during conversion.

```sh
ctest --test-dir build --output-on-failure
# Optional numerical compatibility comparison against the source snapshot:
ALTANN_REFERENCE_BINARY=/path/to/snapshot/uniann ctest --test-dir build --output-on-failure
```

To build a release package locally, configure with `-DALTANN_RELEASE_BUILD=ON`,
install into a staging prefix, then run:

```sh
python3 scripts/package_release.py --prefix /path/to/staging \
  --output dist/altann-linux-amd64.tar.gz
```

The packaging step copies the OpenMP runtime and its notices, and adjusts
library paths on macOS. Keep the extracted package directory intact.

Code layout: `src/model.cpp` defines scoring; `global.cpp` reconstructs the
reference; `local.cpp` performs K-best searches; `output.cpp` converts paths;
`loaders.cpp` validates input tables; `altann/cli.py` handles existing files;
`altann/preprocess.py` prepares PSAURON and site probabilities for each strand;
`altann/strands.py` handles existing strand-tagged scores; `altann/convert.py` selects
candidates and maps GFF coordinates. See
[validation notes](docs/validation.md) for measured compatibility results and
[input formats](docs/format.md) for the detailed file contract.
