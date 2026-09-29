# Changelog

## 0.2.0

- Process one sequence and direction per invocation; use `--reverse` for input
  that was already reverse complemented. Only output coordinates are reflected.
- Require processed emission and transition scores. Remove PSAURON preprocessing
  and the Perl runtime dependency.
- Read a matching UniAnn GFF and DP/BT log by default. Reconstruct exact scores
  from the traceback and inputs, or use `--rerun-viterbi` without a saved log.
- Retain UniAnn's fixed length regulation and seven-state model.
- Replace the multi-direction manifest tutorial with separate chromosome III
  runs using 34.7 MB of prepared S. pombe inputs.

## 0.1.1

- Add S. pombe example inputs and tutorial.
- Clip local flanks before omitted terminal partial models.
- Validate the native UniAnn binary's CDS-only reference GFF.

## 0.1.0

- Standalone seven-state local K-best decoding from existing UniAnn inputs.
- Configurable K (default 10), flanks, and OpenMP threads.
- Both strands, segment mapping, reference-first intron-chain selection, and
  Dmel-compatible score attributes.
- Linux/macOS amd64/arm64 builds and packages with bundled OpenMP runtimes.
- Numerical compatibility tests and chromosome-4 comparisons against the
  existing Dmel implementation.
