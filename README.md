# AltAnn

Alternative gene annotations through local K-best Viterbi decoding.

AltAnn reads existing UniAnn inputs and produces scored alternative gene models
in GFF3. A C++17 core handles dynamic programming and parallel local searches;
a Python standard-library frontend prepares inputs and exports annotations in
original genome coordinates. No UniAnn executable or new export format is
required.

AltAnn builds on **[UniAnn by Aleksey Zimin](https://github.com/alekseyzimin/UniAnn)**
and inherits its **GPLv3** license. Original author credits are retained;
see [AUTHORS.md](AUTHORS.md) and [source provenance](THIRD_PARTY.md).

## Install a release

Release packages target Linux/macOS on amd64/arm64. Download the matching
archive from [Releases](https://github.com/hibiki-kato/AltAnn/releases), extract
it, and run its `bin/altann` launcher. The native decoder and its OpenMP runtime
are included: **users do not need a compiler or a separate OpenMP installation**.

```sh
tar -xzf altann-linux-amd64.tar.gz
altann-linux-amd64/bin/altann decode --input /path/to/uniann-job \
  --strand plus --k 10 --threads auto --output alternatives.gff3
```

Python 3.9 or later is needed for the frontend. Perl is needed only for original
PSAURON input; retained emission tables avoid that dependency. These are
directory packages containing a native executable and a Python frontend, not
single-file executables. Linux release builds use Ubuntu 22.04 (glibc 2.35);
macOS packages are built on macOS 15. Older systems have not been validated.

## Build from source (developers)

Targets: **Linux and macOS, each on amd64 (x86_64) and arm64 (aarch64)**.
Requirements: a C++17 compiler with OpenMP, CMake, and Python 3.9 or later.
Perl is needed only when regenerating emissions from saved PSAURON scores.

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
  --strand plus --k 10 --flank 1000 --threads auto \
  --output alternatives.gff3
```

Each job contains one strand-oriented FASTA record. The directory reader uses
either retained score tables or saved original inputs:

| Input | Accepted files |
| --- | --- |
| Sequence | One `.fa`, `.fasta`, or `.fna` file |
| Retained emissions | `out.ps.txt` |
| Retained motif scores | `out.gt.txt`, `out.ag.txt`, `out.atg.txt`, `out.stop.txt` |
| Original emission input | `psauron_score.csv`, generated with PSAURON `-a` |
| Original motif input | A single `sites*.tsv` |
| Optional execution log | `run.log`, `uniann.log`, or `*.uniann.log` |
| Optional reference check | A single `*.uniann.gff` |

Ambiguous file selection fails with an explanation. Explicit paths override
discovery: `--fasta`, `--emissions`, `--gt`, `--ag`, `--atg`, `--stop`,
`--psauron`, `--sites`, `--log`, and `--best-gff`.

```sh
bin/altann decode --fasta segment.fa --emissions out.ps.txt \
  --gt out.gt.txt --ag out.ag.txt --atg out.atg.txt --stop out.stop.txt \
  --strand plus --threads 8 --output alternatives.gff3
```

Original score conversion reproduces UniAnn's preprocessing. When supplied,
the log supplies the `Multiplier is ..., Factor is ...` values. Without a log,
AltAnn uses the compatibility wrapper's default multiplier and donor-based
factor calculation. Retained score tables are preferable if a run used modified
preprocessing. A rounded DP dump is not used for numerical scoring.

**The global best is recomputed from the same scoring inputs.** A supplied
UniAnn GFF checks that its exon structures are a subset of the reconstructed
reference, allowing UniAnn's short-CDS filtering. GFF alone is insufficient to
recover all omitted genes or the boundary history. AltAnn does not silently
substitute a supplied GFF for the state path. Use matching inputs and the
compatible UniAnn model described below.

## Both strands and multiple chromosomes

Minus-strand inputs must already be reverse complemented, with the matching
PSAURON and site scores. AltAnn maps the resulting annotations back to the
original genome and preserves CDS phase. `--strand minus` does not reverse
the supplied sequence or score arrays. A `.rc.fa` filename or an unambiguous
`plus`/`minus` directory token can supply the orientation automatically.

For multiple jobs, provide a TSV inventory of existing files:

```text
fasta	strand	score_dir
chr1/plus/chr1.fa	plus	chr1/plus
chr1/minus/chr1.rc.fa	minus	chr1/minus
chr2/plus/chr2.fa	plus	chr2/plus
```

```sh
bin/altann decode --manifest inputs.tsv --threads 8 --output genome.gff3
```

Paths are relative to the manifest. Other supported columns are `emissions`,
`gt`, `ag`, `atg`, `stop`, `psauron`, `sites`, `log`, `best_gff`, `seqid`,
`offset`, `sequence_length`, and `segment` (a unique segment identifier).
Segments run sequentially; local loci within each segment run in parallel.
Results are sorted deterministically by sequence ID, coordinates, strand,
and transcript ID.

Headers of the form `chromosome__start-end-total` carry zero-based half-open
segment coordinates. Other headers use the complete sequence by default;
`--seqid`, `--offset`, and `--sequence-length` specify an explicit mapping.
For overlapping segments, midpoint ownership and internal-edge exclusion use
`--segment-overlap 4000000` and `--segment-margin 20000`. Set these to the actual
segmentation settings; they are not inferred from adjacent sequence files.

## Model and selection

The compatibility model has seven coarse states:
`N, E0, E1, E2, I0, I1, I2`. The emission table has five scores per position:
`N, E0, E1, E2, I`; the intron score is shared by all three intron frames.
An optional final nucleotide column is accepted. Positions in score tables
are zero-based; sparse motif files contain `position score` pairs.

AltAnn retains K paths per coarse state (`--k`, default **10**; use `--k 5`
for a smaller candidate set). Length and frame history accompany
each path to enforce UniAnn's transition rules, but **do not expand the state
space**. This preserves the Dmel compatibility search; it is not a guarantee
of the mathematical global top K over all history-dependent paths.

The default search includes each reference gene plus at most 1,000 intergenic
bases on each side, clipped by neighboring genes and sequence ends. Both local
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
| `kbest_reference_score` | Global best replayed over the same local interval |
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
- `.gff3.json`: input paths, parameters, scaling, counts, and elapsed times.

The main GFF is atomically replaced after successful decoding and conversion.
Input files are never modified.

## Performance and development

`--threads auto` uses CPUs available to the process. `--threads N` provides
explicit control. Memory grows with sequence length and with the number and
length of loci decoded concurrently; use a smaller count on memory-constrained
machines. Auto currently selects CPUs, not a memory budget.

The global pass keeps rolling scores and compact traceback. The local pass
parallelizes independent loci with OpenMP and writes results in fixed order.
The Python frontend uses disk-backed sorting across segments. A segment's
candidate annotations are held in memory during conversion.

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
