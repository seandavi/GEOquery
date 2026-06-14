<!--
Thanks for contributing to GEOquery! Please complete the checklist below.
See CONTRIBUTING.md for the full development workflow and conventions.
-->

## What does this PR do?

<!-- One or two sentences. Link the issue it closes, e.g. "Fixes #123". -->

Fixes #

## Checklist

- [ ] **NEWS** — added a user-facing entry to `NEWS.md` under `# GEOquery (development version)`, ending with the PR link, e.g. `* `getGEO()` now ... (#NNN).` (Skip only for CI/internal-only changes.)
- [ ] **Tests** — added or updated tests covering the change. Network-dependent tests are gated (`skip_if_offline()` / fixtures) where applicable.
- [ ] **Docs** — ran `devtools::document()`; `man/` and `NAMESPACE` are regenerated and committed (do not hand-edit).
- [ ] **Version** — bumped `Version:` (`z`) and `Date:` in `DESCRIPTION` (Bioconductor convention; skip for build-excluded changes).
- [ ] **Checks** — `R CMD check` and `BiocCheck::BiocCheck()` pass locally (or the required CI legs are green).
- [ ] **ADR** — for a non-trivial architectural decision, added an `adr/NNNN-*.md` (see `adr/template.md`).

## Notes for reviewers

<!-- Anything reviewers should focus on, trade-offs, or follow-ups. -->
