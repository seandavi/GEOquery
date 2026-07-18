# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

GEOquery is a Bioconductor R package that downloads and parses data from
NCBI Gene Expression Omnibus (GEO) into Bioconductor data structures. It
is the bridge between GEO’s SOFT/Series-Matrix file formats and
Bioconductor classes (`ExpressionSet`, `SummarizedExperiment`).

## Commands

R package — no build step for editing; reload and test via
`devtools`/`testthat`.

``` r

devtools::load_all(".")        # reload package after edits
devtools::document()           # regenerate man/ and NAMESPACE from roxygen2 (Config/roxygen2/version: 8.0.0)
devtools::test()               # run full testthat suite
devtools::test(filter = "GSE") # run one test file (tests/testthat/test_GSE.R)
testthat::test_file("tests/testthat/test_GSE.R")  # run a single file directly
```

If `devtools` is unavailable, `pkgload::load_all(".")`,
`roxygen2::roxygenize(".")`, and `testthat::test_dir("tests/testthat")`
cover load/document/test without it.

### Local R environment setup (from a bare R)

On a fresh machine R may have **no** GEOquery dependencies (only
base/recommended packages). Setup that works here:

- **Library location:** use R’s *default* user library
  (`~/R/<platform>-library/<ver>`). Do **not** set `R_LIBS_USER` to a
  custom path — just create the default dir if missing
  (`dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE)`); R auto-adds
  it to [`.libPaths()`](https://rdrr.io/r/base/libPaths.html) for every
  session. The system lib (`/usr/lib64/R/library`) is not writable.
- **System dev libraries** are needed before some R packages build (this
  is RHEL/el9; passwordless `sudo dnf` is available): `libxml2-devel`
  (for `xml2`, a hard dep), plus
  `fontconfig-devel freetype-devel harfbuzz-devel fribidi-devel libpng-devel`
  (for the `ragg`/`textshaping` doc-tooling chain), and
  `libcurl-devel openssl-devel`. `pak` will also auto-run `dnf` for a
  package’s declared SystemRequirements.
- **Installer:** prefer `pak::pkg_install(...)` — fast, parallel, robust
  locking, resolves Bioconductor automatically, and installs
  SystemRequirements. Bootstrap it with `install.packages("pak")`.
  `BiocManager::install()` with `Ncpus>1` was prone to `00LOCK-*` lock
  cascades here; if you hit one, `rm -rf <lib>/00LOCK-*` and retry (pak
  avoids this).
- **Skip the doc-website tooling** (`devtools` → `pkgdown`/`ragg`)
  unless building the site: `pkgload` (load_all), `roxygen2` (document),
  and `testthat` (test) are enough and avoid the font stack.
- **Match the roxygen version:** docs are generated with the version in
  `DESCRIPTION`’s `Config/roxygen2/version` (currently 8.0.0). A
  mismatched roxygen2 reformats every `man/*.Rd`, creating churn —
  install the pinned version (`pak::pkg_install("roxygen2@8.0.0")`)
  before `roxygenize()`. **Never hand-edit `man/*.Rd` or `NAMESPACE`** —
  they are generated.
- **Running integration tests locally:** the network tests are gated by
  `skip_if_no_integration()` *and* `skip_on_cran()`. Set both
  `NOT_CRAN=true` and `GEOQUERY_INTEGRATION=true` to actually run them
  (devtools::test sets `NOT_CRAN` for you; raw `testthat::test_*` does
  not).

Full checks (match CI in `.github/workflows/R-CMD-check.yaml`):

``` sh
R CMD build .
R CMD check GEOquery_*.tar.gz
Rscript -e 'BiocCheck::BiocCheck()'   # Bioconductor-specific checks
Rscript -e 'lintr::lint_package()'    # config in .lintr
```

Docker dev environment: `Dockerfile` builds on
`bioconductor/bioconductor_docker:devel`.

### Important: tests hit the live network

The testthat suite makes real HTTP calls to NCBI GEO FTP/CGI endpoints —
there are no recorded fixtures. Tests download specific accessions
(e.g. GSE/GSM/GPL/GDS IDs) and assert on parsed structure. Failures may
reflect NCBI outages or GEO-side format changes, not code regressions.
When a parsing bug is reported it is usually tied to a specific
accession; reproduce with `getGEO("GSExxxxx")` against the real record.

## Architecture

### Entry point and dispatch

[`getGEO()`](http://seandavi.github.io/GEOquery/reference/getGEO.md)
(`R/getGEO.R`) is the sole user entry point. Flow:

1.  If `GEO` is a GSE **and** `GSEMatrix=TRUE` (default) →
    `getAndParseGSEMatrices()` → fast Series Matrix path → list of
    `ExpressionSet`s. This bypasses SOFT parsing entirely.
2.  Otherwise →
    [`getGEOfile()`](http://seandavi.github.io/GEOquery/reference/getGEOfile.md)
    downloads the SOFT file →
    [`parseGEO()`](http://seandavi.github.io/GEOquery/reference/parseGEO.md)
    dispatches on entity type.

[`parseGEO()`](http://seandavi.github.io/GEOquery/reference/parseGEO.md)
(`R/parseGEO.R`) reads the first entity marker (`findFirstEntity`) and
switches to `parseGSM` / `parseGSE` / `parseGDS` / `parseGPL`, or to
`parseGSEMatrix` when given a Series Matrix file (entity `0`).

There are **two distinct parse paths** that produce different return
types — keep them separate when fixing bugs: - **SOFT format** →
GEOquery S4 objects (`GSE`/`GSM`/`GPL`/`GDS`). - **Series Matrix**
(`parseGSEMatrix`) → Bioconductor `ExpressionSet`. This is the default
for GSEs and is orders of magnitude faster.

### S4 class hierarchy (`R/classes.R`)

- `GEOData` — base, holds `header` (metadata list). Accessors:
  [`Meta()`](http://seandavi.github.io/GEOquery/reference/GEOData-accessors.md),
  [`Accession()`](http://seandavi.github.io/GEOquery/reference/GEOData-accessors.md).
- `GEODataTable` — `columns` (descriptions) + `table` (data)
  data.frames. Accessors:
  [`Table()`](http://seandavi.github.io/GEOquery/reference/GEOData-accessors.md),
  [`Columns()`](http://seandavi.github.io/GEOquery/reference/GEOData-accessors.md).
- `GSM`, `GPL`, `GDS` — extend `GEOData`, each carry a `dataTable` slot.
- `GSE` — container only: `header` + `gsms` (list of `GSM`) + `gpls`
  (list of `GPL`). Accessors:
  [`GSMList()`](http://seandavi.github.io/GEOquery/reference/GEOData-accessors.md),
  [`GPLList()`](http://seandavi.github.io/GEOquery/reference/GEOData-accessors.md).

When adding accessors: define a `setGeneric` here and `#' @export` the
`setMethod` so roxygen exports it.

### Parsing internals (`R/parseGEO.R`, the largest file)

- Files are read with `data.table::fread(sep="")` into a single
  character vector, then dissected by regex. Note
  `na.strings`/`na = .na_strings` handling and `fill=TRUE` are
  load-bearing for malformed GEO files (see NEWS — GSE425 fix extracts
  `!Sample_` lines by pattern match rather than positional reads).
- `parseGeoMeta()` / `parseGeoColumns()` split `!key = value` and
  `#col = desc` lines.
- `.genericGEOTableParser()` locates `!..._table_begin`/`_table_end`
  markers to separate metadata from the data table.
- `fastTabRead()` guesses column classes from the first ~100 rows for
  speed — a deliberate correctness/speed tradeoff; representative-sample
  assumptions can bite on heterogeneous tables.
- [`parseGSEMatrix()`](http://seandavi.github.io/GEOquery/reference/parseGSEMatrix.md)
  builds the `ExpressionSet`: header → `MIAME`, `!Sample_` lines →
  phenoData, `characteristics_ch1/ch2` key:value pairs are unpacked into
  pData columns (`parseCharacteristics`), and the GPL is fetched for
  featureData unless `getGPL=FALSE`.

### Other subsystems

- `R/getGEOfile.R` — URL construction (accession → FTP path via `nnn`
  stubbing) and `downloadFile()` (httr2-based).
- `R/getGEOSuppFiles.R` — supplementary file listing/download;
  [`getDirListing()`](http://seandavi.github.io/GEOquery/reference/getDirListing.md)
  scrapes FTP directory HTML.
- `R/rnaseq.R` — NCBI-computed RNA-seq quantification support
  (`getRNASeqData`, `hasRNASeqQuantifications`); scrapes GEO download
  pages with `rvest`/`xml2` to find raw-count and annotation URLs.
- `R/searchGEO.R` — Entrez search via `rentrez` (`searchGEO`,
  `searchFieldsGEO`).
- `R/GDS2MA.R` — convert GDS to limma `MAList` / `ExpressionSet`.

### Runtime options (`R/zzz.R`)

Set on load: `download.file.method.GEOquery = "auto"` and
`GEOquery.inmemory.gpl = FALSE`. The latter gates an in-memory GPL cache
(`GPLcache` env in `parseGEO.R`) keyed by file mtime — off by default.

## Conventions

- Object naming: `snake_case`, `camelCase`, or symbol allowed (`.lintr`
  sets `object_name_linter`). Existing code mixes `camelCase` functions
  (`parseGSEMatrix`) with `.dotPrefixed` internals (`.parseGSMTxt`,
  `.read_lines`).
- Roxygen2 with markdown (`Roxygen: list(markdown = TRUE)`); `man/` and
  `NAMESPACE` are generated — never hand-edit, run
  `devtools::document()`.
- Add user-facing changes to `NEWS.md` under the current unreleased
  version heading.
- Vignettes are Quarto (`.qmd`, `VignetteBuilder: quarto`), not Rmd, and
  are **precompiled**: edit `vignettes/<name>.qmd.orig` (live code),
  then run `Rscript dev/precompute-vignettes.R` to regenerate the
  shipped static `vignettes/<name>.qmd` (real GEO output baked in, no
  network at build). Never hand-edit the generated `.qmd`; commit both
  files. See CONTRIBUTING.md.
- Branching follows Bioconductor: `devel` is the main working branch.

## Architecture Decision Records

Significant architectural decisions are recorded as ADRs in `adr/`
(`NNNN-title.md`, numbered sequentially, following `adr/template.md`).
See `adr/0001-record-architecture-decisions.md` for the rationale and
rules. The directory is excluded from the package build via
`.Rbuildignore`.

When making a non-trivial architectural choice (a new parse path, a
changed return type, a dependency swap, a speed/correctness tradeoff),
write an ADR:

1.  Copy `adr/template.md` to the next number.
2.  Fill in Context / Decision / Consequences / Alternatives.
3.  ADRs are immutable once `accepted` — to revise, write a new ADR that
    supersedes it and set the old one’s status to
    `superseded by ADR-XXXX`.
