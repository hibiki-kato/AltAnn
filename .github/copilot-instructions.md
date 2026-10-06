# AltAnn review instructions

Review changes for correctness, reproducibility, and readable English code.
Prioritize actionable defects over style-only remarks.

Write clear, natural English. Do not use en dashes, em dashes, canned AI
phrasing, promotional language, or unnecessary contrast constructions in prose.
Use comments to explain purpose, assumptions, and boundary cases. Keep existing
upstream license text unchanged.

- The model has exactly seven coarse states. Do not introduce metadata-expanded
  state buckets or change numerical transition rules without an explicit design
  change and compatibility evidence.
- K defaults to 10 and is configurable. It counts retained paths before
  annotation deduplication, not a fixed number of exported isoforms.
- Internal coordinates are zero-based half-open. GFF coordinates are one-based
  inclusive. In both-strand mode, site probabilities use one-based original
  segment coordinates and PSAURON supplies six frame arrays. The frontend
  prepares reverse inputs and processed scores. Legacy `--reverse`
  inputs are already reverse complemented. Verify reflection,
  segment offsets, CDS phase, and ownership at overlapping segment boundaries.
- Reference rank 0 is the actual reconstructed UniAnn best. Alternative deltas
  compare path scores over the same local interval with identical initial
  history. Scores are not probabilities. Preserve full serialization precision.
- Keep one candidate per intron chain within a locus, with reference precedence.
  Do not deduplicate across unrelated loci or silently replace filtered references.
- Preserve deterministic output across thread counts. Exceptions must not escape
  OpenMP regions. Review bounds, allocation growth, and malformed score inputs.
- Never execute input text as shell commands or overwrite input data. Keep
  subprocess arguments as lists. Public workflows must not execute untrusted PR
  code with write tokens or secrets.
- Support Linux/macOS on amd64/arm64. Release packages must use their bundled
  OpenMP runtime, without requiring a user's compiler or Homebrew installation.
- Preserve UniAnn's GPLv3 license and Aleksey Zimin's authorship credits.
- Require meaningful regression coverage for numerical, coordinate, parallel,
  input-format, and packaging changes. Explain biological and compatibility
  assumptions in comments and docs, all in English.
