# Changelog

## Unreleased

- Apply `--short-cds-filter` before deduplicating, so a removed model no longer
  hides a later model with the same intron chain. Treat a single-exon model as a
  duplicate only when it overlaps a retained single-exon model of the locus, not
  whenever the locus already has one.

- Add `altann-core --uniann-gff`, which decodes the global best path and prints it
  as the UniAnn binary does, without UniAnn's memory use or DP/BT log.

- Trace all K paths back in one pass, store them run-length encoded, and spill
  tracebacks that do not fit the shared `--traceback-memory` budget to a temporary
  file, so long windows at high K fit in memory at any thread count. Output is
  unchanged.

- Read emission scores on demand from a position index instead of holding them
  in memory, keep motif scores sparse, and store the global path and N-prefix
  metadata in compact form. Peak memory on Dmel chrX falls from 2.17 GB to
  0.37 GB; output is unchanged. Emission inputs must now be seekable files.

- Add `altann decode --raw` to export every gene of every decoded path, without
  the CDS-overlap filter or annotation and intron-chain deduplication. Default
  output is unchanged.

- Add `altann score` to replay supplied GFF3/GTF CDS models under the native
  model and report their delta to the global-best path on the decoder window,
  or the first forbidden transition.

- Report each candidate's donor/acceptor `splice_site_jaccard` against its
  global-best reference in the diagnostic TSV and as the GFF3 attribute
  `kbest_splice_jaccard`.

- Retain local alternatives only when their CDS overlaps their own global-best
  reference CDS; exclude independent flank and intronic ORFs before ranking
  without changing the reference models, path scores or Viterbi decoding.
- Complete real CNN-Mamba → UniAnn → AltAnn validation on all S. pombe nuclear
  chromosomes and both strands; add evaluation and spliced-CDS audit scripts
  with real-DNA regressions for shared upstream biological limitations.
- Validate the local UniAnn header/frame/calibration fixes with corrected
  score fixtures and a real S. pombe reference-coordinate audit; document the
  remaining GPU/model-output blocker separately from prediction accuracy.
- Consume UniAnn `6fe3d46` coding-frame emissions through the existing processed
  score interface, without a model runtime or additional calibration.
- Prefer processed scores during directory discovery even with `-a` and a
  leftover PSAURON CSV; retain explicit legacy raw-input support.
- Document local-to-absolute frame conversion, upstream score calibration,
  and the current oriented UniAnn workflow.
- Add upstream-generated coding-frame fixtures and regressions for all frames,
  splice phases, score replay, alternative acceptors, and reverse coordinates.
- Inspect the actual cnn-mamba exporter and calibration, add reproducible
  downstream contract diagnostics, and document upstream header handling,
  one-frame indexing offset, and repeated temperature transformations.

## 0.3.0

- Decode both strands in one invocation with `-a` / `--all-prob`, using
  the original FASTA, six-frame PSAURON CSV, and both-strand site probabilities.
- Prepare emission and transition scores internally with UniAnn's rules,
  using Python's standard library. Accept `-f`, `-p`, `-s`, and `-m` inputs.
- Recompute both baseline paths independently and export one GFF3 and diagnostics
  table. Preserve the seven-state decoder, score precision, and CDS phase.
- Retain processed UniAnn GFF/traceback inputs and the existing `--reverse` mode.
- Update the chromosome III tutorial to produce one combined annotation.
- Fix Makefile installation after removal of the Perl preprocessor.

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
