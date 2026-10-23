# Source provenance

AltAnn inherits UniAnn's GNU GPL version 3 license. [LICENSE](LICENSE) is an
unchanged copy of the upstream license, not a newly selected license.

The seven-state scoring model, transcript conversion, and coarse-state local
K-best search are adapted from **[UniAnn by Aleksey Zimin](https://github.com/alekseyzimin/UniAnn)**,
through [Hibiki Kato's UniAnn fork](https://github.com/hibiki-kato/UniAnn).
The compatibility baseline is the UniAnn source used in the Dmel
`Dmel_k20_flank1000_chain` experiment: commit
`d0722611ab47f0f59d1391fb015b0248194ca836`, with the experiment's bounded-flank,
inherited-boundary-history, and start-codon boundary patches. AltAnn separates
these routines into modules, uses a rolling global score array, and parallelizes
independent local searches. It does not use metadata-expanded states.

GFF export is adapted from the same Dmel experiment's conversion script.
Credits for the original UniAnn implementation remain with its authors and
contributors. The probability preprocessing in `altann/preprocess.py` is adapted
from UniAnn's `scripts/preprocess_psauron_scores.pl` and `scripts/uniann.sh` at
commit `91477a69e1a949fed082c4662b348d9a7a91ca2b`. It preserves the emission
transformation, stop-codon handling, rounding, and site probability scaling.
The Python implementation requires no Perl runtime or UniAnn executable.

Processed-score compatibility was checked against upstream UniAnn commit
`6fe3d46efacf1f900789a4cfa6d103c69dfe6ac0`. The native `src/uniann.cpp`
is unchanged from `91477a6`; its new `scripts/convert_coding_frame.pl` maps
local coding-frame probabilities to the same emission columns. That converter's
temporary calibration was not copied into AltAnn; local UniAnn fixes now
remove it and correct header handling and local-frame indexing.
New regression inputs are produced by the upstream converter and decoder,
as documented in `tests/data/README.md`.

The predictor contract was also reviewed against Bowen's
[cnn-mamba](https://github.com/bowencccc/cnn-mamba) at
`809f74379efccf70ee7fd876dc95aabe5a4b0ea7`. The `bowen-contract` fixtures
exercise its actual phase-label, temperature, and TSV export code on synthetic
logits. None of that model or calibration code is copied into AltAnn's runtime.
The original model/converter discrepancies and local UniAnn fixes are
documented in `docs/format.md`. The `bowen-corrected` fixtures replay the same
export through corrected source scripts; these are integration regressions,
not trained-model accuracy measurements. Real S. pombe validation status is
recorded in `docs/spom-validation.md`.

## Credits

- **Aleksey Zimin and UniAnn contributors:** original UniAnn gene prediction
  implementation, scoring rules, and preprocessing code.
- **Hibiki Kato:** local K-best extensions and Dmel evaluation work, and the
  AltAnn standalone extraction, interfaces, packaging, and maintenance.

This attribution does not transfer ownership of upstream work. Existing
upstream notices are preserved. See the original repository history for the
full contributor record. The GPL document's Free Software Foundation copyright
notice applies to the license text and is not a software authorship credit.

AltAnn requires no UniAnn executable, neural-network runtime, PSAURON executable,
gffread, or gffcompare when decoding saved inputs. The S. pombe genome and
precomputed scores in `example/data/` are documented in
[the example tutorial](example/README.md). Generated annotation results are
excluded from version control.
