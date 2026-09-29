# Changelog

## 0.1.0

- Standalone seven-state local K-best decoding from existing UniAnn inputs.
- Configurable K (default 10), flanks, and OpenMP threads.
- Both strands, segment mapping, reference-first intron-chain selection, and
  Dmel-compatible score attributes.
- Linux/macOS amd64/arm64 builds and packages with bundled OpenMP runtimes.
- Numerical compatibility tests and chromosome-4 comparisons against the
  existing Dmel implementation.
