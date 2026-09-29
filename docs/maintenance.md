# Maintenance and releases

All source changes go through pull requests to `main`. The repository ruleset
requires four native build/test/package checks, an up-to-date branch, and
resolved review conversations. Force pushes and branch deletion are blocked.
Copilot reviews new PRs and subsequent pushes, with project-specific guidance
in `.github/copilot-instructions.md`. Human approval count is zero so a solo
maintainer is not blocked from merging their own reviewed PR; passing CI and
conversation resolution remain required.

The GitHub Actions token defaults to read-only. Only the release job receives
contents write permission. Pull-request jobs have no release permission.
Dependency updates for Actions are opened monthly by Dependabot.

## Versioned releases

1. Update `altann/__init__.py`, `CMakeLists.txt`, and `CHANGELOG.md` in a PR.
2. Merge after all four platform jobs pass.
3. The main-branch pipeline rebuilds and packages all four platforms.
4. If that version has no tag, the release job creates `vX.Y.Z` at the tested
   commit and publishes the matching archives and checksums.

Existing tags/releases are never overwritten. Commits that do not change the
version still run all checks but do not create duplicate releases. Tagging and
publishing happen in the same workflow, avoiding reliance on workflow-trigger
events created by `GITHUB_TOKEN`.

Release packages contain the native decoder, Python frontend, compatibility
preprocessor, runtime notices, and bundled OpenMP library. End users need
Python 3.9+ (and Perl only for raw PSAURON input), not a compiler or OpenMP
development installation. Source code remains available with each release tag.

## Package managers

GitHub Releases is the initial distribution channel. A separately maintained
Homebrew tap can later point at these releases; acceptance into Homebrew's
official repository is a different process. Bioconda can be considered once
the interface and release process are stable. Neither registration is required
for users to run the release packages.
