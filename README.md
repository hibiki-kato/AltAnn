# AltAnn

Alternative gene annotations through local K-best Viterbi decoding.

AltAnn reads existing UniAnn inputs and produces scored alternative gene models
in GFF3. A C++17 core handles dynamic programming and parallel local searches;
a Python standard-library frontend reads processed score tables and exports
annotations in original genome coordinates. Each invocation decodes one FASTA
record in one orientation. No UniAnn executable is required.

AltAnn builds on **[UniAnn by Aleksey Zimin](https://github.com/alekseyzimin/UniAnn)**
and inherits its **GPLv3** license. Original author credits are retained;
see [AUTHORS.md](AUTHORS.md) and [source provenance](THIRD_PARTY.md).

Try the [S. pombe example](example/README.md) to decode chromosome III using
supplied score files. Forward and reverse inputs are separate runs.

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
PSAURON CSV and site probability tables must be processed by UniAnn beforehand;
AltAnn consumes the resulting emission and transition scores.

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

## Reverse inputs

**Before using `--reverse`, you must supply a reverse-complemented FASTA.**
Reverse complement means reversing the sequence and replacing each base with
its complement (A to T, T to A, C to G, G to C). Reversing the base order alone is insufficient.
AltAnn does not perform this conversion for you.

One invocation processes one orientation. To annotate the opposite strand,
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

Run additional sequences or orientations separately. Merge their GFF files as
a separate downstream step if a combined annotation is needed.

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

Reference transcripts are retained with rank 0. Candidates use the original
local decoder ranks, up to K. Intron chains are deduplicated within each locus;
the reference takes precedence over alternatives sharing its chain. Empty
chains (single-exon models) also deduplicate within the locus. K is the number
of retained local paths, not a promise of K exported isoforms. A path may
contain several genes. Selection happens before optional output filtering.

`--short-cds-filter` reproduces the Dmel export rule: omit transcripts with at
most two exons and at most 200 coding bases. It is disabled by default so this
filter is an explicit choice.

The model follows the coarse-state UniAnn snapshot documented in
[THIRD_PARTY.md](THIRD_PARTY.md). Other UniAnn versions or customized scoring
rules may produce a different global best and must be checked before use.

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

## Performance and development

`--threads auto` uses CPUs available to the process. `--threads N` provides
explicit control. Memory grows with sequence length and with the number and
length of loci decoded concurrently; use a smaller count on memory-constrained
machines. Auto currently selects CPUs, not a memory budget.

The global pass keeps rolling scores and compact traceback. The local pass
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
`altann/convert.py` selects candidates and maps GFF coordinates. See
[validation notes](docs/validation.md) for measured compatibility results and
[input formats](docs/format.md) for the detailed file contract.
