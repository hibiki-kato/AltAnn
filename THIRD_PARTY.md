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

`vendor/preprocess_psauron_scores.pl` is the UniAnn compatibility preprocessor.
It is retained to reproduce existing PSAURON-to-emission conversion and requires
Perl only when original PSAURON input is supplied. Its numerical transformations
are not rewritten in Python. GFF export is adapted from the same Dmel
experiment's conversion script. Credits for the original UniAnn implementation
remain with its authors and contributors.

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
gffread, or gffcompare when decoding saved inputs. External sequence data and
research result files are not included in this repository.
