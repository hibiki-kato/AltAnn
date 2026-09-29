# Validation

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

## Timing

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
