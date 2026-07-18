# Contributing to GEOquery

Thanks for contributing. This guide covers the development workflow, the
changelog convention, and the checks a pull request must pass. GEOquery is a
[Bioconductor](https://bioconductor.org) package; `devel` is the working branch.

## Issues vs. pull requests

Not every change needs an issue — file one where it helps users and
contributors, not as ceremony.

- **Bug, feature, or behavior change → open an issue first** (or link an
  existing one), then reference it from the PR with `Fixes #123`. These are
  user-facing, so they also get a `NEWS.md` bullet that links both the issue
  and the PR. For bugs, include a **GEO accession + reproducible example +
  `sessionInfo()`** (the bug issue form prompts for these).
- **Chore, CI, docs, or refactor → a pull request alone is fine.** No issue
  needed; use a `chore:` / `ci:` / `docs:` / `refactor:` commit prefix.

Planned work is tracked on the [`3.0` milestone](https://github.com/seandavi/GEOquery/milestones);
see `ROADMAP.md` for the broader plan.

## Development workflow

1. Start work off `devel`. With **jujutsu** (`jj`, the repo is colocated so git
   still works): `jj new devel`, make the change, then `jj describe -m "..."`.
   With plain git: `git switch -c fix/issue-123`. Either way `devel` is `trunk()`.
2. Make the change. Keep PRs focused — one logical change per PR.
3. Reload and test locally:
   ```r
   devtools::load_all(".")
   devtools::test()                 # or devtools::test(filter = "GSE") for one file
   ```
4. Regenerate docs if you touched roxygen:
   ```r
   devtools::document()             # updates man/ and NAMESPACE — never hand-edit these
   ```
5. Update `NEWS.md` (see below) and bump the version (see below).
6. Push and open a pull request against `devel`. With `jj`: `jj git push --named fix/issue-123=@`
   (creates, tracks, and pushes the bookmark). With git: `git push -u origin fix/issue-123`.
   The PR template's checklist is the contract; CI runs `R CMD check` across a platform matrix.
7. Merge when the required checks are green. (`devel` is branch-protected; the
   required legs are the Ubuntu and macOS `R CMD check` jobs.)

## jj for git users (optional)

**Entirely optional.** The repo is a colocated [jujutsu](https://jj-vcs.github.io/jj/)
+ git checkout, so every `git` command still works and the remote only ever sees
plain git branches — use jj only if you want to. To set it up: `jj git init --colocate`
in your clone. Then the common steps map like this (`devel` is `trunk()`, jj
bookmarks are git branches):

| git | jj |
|-----|-----|
| `git switch -c fix/x devel` | `jj new devel` |
| `git status` / `git log --oneline` | `jj st` / `jj log` |
| `git commit -m "..."` | `jj describe -m "..."` (the working copy *is* the commit — no staging) |
| `git rebase devel` | `jj rebase -o devel` |
| `git push -u origin fix/x` | `jj git push --named fix/x=@` |
| `git push` (updates) | `jj git push` |

Push still goes through a bookmark → PR against `devel`; nothing about the
Bioconductor side changes.

Safety net: `jj op log` shows every operation jj has performed, and
`jj op restore <id>` (or `jj undo` for the last one) rewinds the *whole repo* to
that state — so a botched rebase or an accidental `jj abandon` is one command to
recover.

## Vignettes are precompiled

The vignettes hit the live NCBI GEO network, which must **not** happen during
`R CMD build`/`check` on CRAN/Bioconductor build machines. So each vignette is
authored in a `vignettes/<name>.qmd.orig` source (with live `eval: true` code)
and **precompiled** into the shipped static `vignettes/<name>.qmd`, which has the
real output baked in and executes nothing at build time. (This is the knitr
`*.Rmd.orig` pattern; Quarto's `freeze` does *not* help here, because the vignette
engine renders each file individually and freeze only applies to full-project
renders.)

Do **not** edit `vignettes/*.qmd` directly — they are generated. To change a
vignette:

1. Edit the `vignettes/<name>.qmd.orig` source.
2. Regenerate (with network access and the `Suggests` deps installed):
   ```sh
   Rscript dev/precompute-vignettes.R            # all vignettes
   Rscript dev/precompute-vignettes.R rnaseq     # just one
   ```
3. Commit **both** the `.qmd.orig` source and the regenerated `.qmd`.

Refresh whenever you change a vignette's code, or periodically to pick up
GEO-side changes. Chunks that can't run at precompile time (private-token
examples, very large single-cell downloads, Seurat coercions, pure pseudo-code)
are marked `#| eval: false` in the source and carry hand-written output.

## NEWS / changelog convention

We follow the [tidyverse NEWS style](https://style.tidyverse.org/news.html).
Every user-facing change gets one bullet in `NEWS.md` under the top
`# GEOquery (development version)` header (this header is renamed to the release
version at release time).

Rules:

- Write a **complete sentence** describing the change from the **user's**
  perspective. Put function/argument names in backticks.
- **End every bullet with the PR link**, and credit the contributor:
  ```markdown
  * `getGEO(parseCharacteristics = FALSE)` no longer parses characteristics;
    the flag is now threaded through `parseGSEMatrix()` (#166, @reporter).
  ```
- Group bullets under `## Breaking changes`, `## New features`, or
  `## Bug fixes` only once there are enough to warrant it; otherwise a flat list
  is fine. Most recent entries at the top.
- The `#NNN` and `@user` references autolink in the pkgdown changelog and on
  GitHub (the repo URL comes from `DESCRIPTION`).
- **Skip NEWS** for CI-only or internal-only changes (nothing user-visible).

## Versioning

Bump `Version:` and `Date:` in `DESCRIPTION` on every code PR, following the
Bioconductor convention (increment the `z` in `x.y.z` during the `devel` cycle).
Build-excluded changes (anything matched by `.Rbuildignore`, e.g. `.github/`,
`CONTRIBUTING.md`, `ROADMAP.md`) do not require a version bump.

## Checks before merge

- `R CMD check` — no errors or warnings.
- `BiocCheck::BiocCheck()` — Bioconductor-specific requirements.
- `lintr::lint_package()` — style (config in `.lintr`).

## Architecture Decision Records

Non-trivial architectural decisions are recorded as ADRs in `adr/` (see
`adr/template.md` and `CLAUDE.md`). If your PR changes a parse path, a return
type, a dependency, or a notable trade-off, add a numbered ADR.

## Reporting bugs

Open an issue with a **GEO accession**, a **minimal reproducible example**, and
`sessionInfo()`. Most GEOquery bugs are specific to how one record is formatted
on NCBI, so the accession is essential for reproduction.
