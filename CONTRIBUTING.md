# Contributing to GEOquery

Thanks for contributing. This guide covers the development workflow, the
changelog convention, and the checks a pull request must pass. GEOquery
is a [Bioconductor](https://bioconductor.org) package; `devel` is the
working branch.

## Issues vs. pull requests

Not every change needs an issue — file one where it helps users and
contributors, not as ceremony.

- **Bug, feature, or behavior change → open an issue first** (or link an
  existing one), then reference it from the PR with `Fixes #123`. These
  are user-facing, so they also get a `NEWS.md` bullet that links both
  the issue and the PR. For bugs, include a **GEO accession +
  reproducible example +
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)** (the bug
  issue form prompts for these).
- **Chore, CI, docs, or refactor → a pull request alone is fine.** No
  issue needed; use a `chore:` / `ci:` / `docs:` / `refactor:` commit
  prefix.

Planned work is tracked on the [`3.0`
milestone](https://github.com/seandavi/GEOquery/milestones); see
`ROADMAP.md` for the broader plan.

## Development workflow

1.  Branch off `devel` (e.g. `git switch -c fix/issue-123`).

2.  Make the change. Keep PRs focused — one logical change per PR.

3.  Reload and test locally:

    ``` r

    devtools::load_all(".")
    devtools::test()                 # or devtools::test(filter = "GSE") for one file
    ```

4.  Regenerate docs if you touched roxygen:

    ``` r

    devtools::document()             # updates man/ and NAMESPACE — never hand-edit these
    ```

5.  Update `NEWS.md` (see below) and bump the version (see below).

6.  Open a pull request against `devel`. The PR template’s checklist is
    the contract; CI runs `R CMD check` across a platform matrix.

7.  Merge when the required checks are green. (`devel` is
    branch-protected; the required legs are the Ubuntu and macOS
    `R CMD check` jobs.)

## NEWS / changelog convention

We follow the [tidyverse NEWS
style](https://style.tidyverse.org/news.html). Every user-facing change
gets one bullet in `NEWS.md` under the top
`# GEOquery (development version)` header (this header is renamed to the
release version at release time).

Rules:

- Write a **complete sentence** describing the change from the
  **user’s** perspective. Put function/argument names in backticks.

- **End every bullet with the PR link**, and credit the contributor:

  ``` markdown
  * `getGEO(parseCharacteristics = FALSE)` no longer parses characteristics;
    the flag is now threaded through `parseGSEMatrix()` (#166, @reporter).
  ```

- Group bullets under `## Breaking changes`, `## New features`, or
  `## Bug fixes` only once there are enough to warrant it; otherwise a
  flat list is fine. Most recent entries at the top.

- The `#NNN` and `@user` references autolink in the pkgdown changelog
  and on GitHub (the repo URL comes from `DESCRIPTION`).

- **Skip NEWS** for CI-only or internal-only changes (nothing
  user-visible).

## Versioning

Bump `Version:` and `Date:` in `DESCRIPTION` on every code PR, following
the Bioconductor convention (increment the `z` in `x.y.z` during the
`devel` cycle). Build-excluded changes (anything matched by
`.Rbuildignore`, e.g. `.github/`, `CONTRIBUTING.md`, `ROADMAP.md`) do
not require a version bump.

## Checks before merge

- `R CMD check` — no errors or warnings.
- `BiocCheck::BiocCheck()` — Bioconductor-specific requirements.
- `lintr::lint_package()` — style (config in `.lintr`).

## Architecture Decision Records

Non-trivial architectural decisions are recorded as ADRs in `adr/` (see
`adr/template.md` and `CLAUDE.md`). If your PR changes a parse path, a
return type, a dependency, or a notable trade-off, add a numbered ADR.

## Reporting bugs

Open an issue with a **GEO accession**, a **minimal reproducible
example**, and
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). Most
GEOquery bugs are specific to how one record is formatted on NCBI, so
the accession is essential for reproduction.
