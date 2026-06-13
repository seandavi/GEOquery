# GEOquery Improvement Roadmap

This roadmap synthesizes six parallel audits (bugs, features, UX, docs, CI/CD, triage automation) into a single, prioritized plan for the GEOquery package. It leads with the highest-leverage quick wins, groups the remaining work by domain, gives special depth to the issue-triage/PR-automation design (the maintainer's focus area), and ends with a phased sequencing plan. The throughline: a handful of small forwarding/validation fixes unblock multiple open issues, while two structural investments — **offline test fixtures** and a **layered triage/automation system** — are the foundations everything else depends on.

---

## Release plan: GEOquery 3.0 & paper

This roadmap doubles as the release plan. The package is on the `2.77.x` Bioconductor devel line (the authoritative `DESCRIPTION` version), so the next release is `2.78.0` by default. **"GEOquery 3.0" is a project/paper codename, not the current version number.** To release as `3.0.0`, the devel version is set to `2.99.z` (Bioconductor's pre-3.0 convention) **deliberately, near the actual release**, once the headline features have shipped — not before. The goal is to make 3.0 a **milestone branded on *shipped* code**, not an aspirational version bump.

**3.0 theme — modern object model + modern data types.** SummarizedExperiment/SingleCellExperiment return types and single-cell support are the headline; bug fixes, CI, and triage automation are supporting cast.

### 3.0 batch (must ship before cutting 3.0)

The minimum coherent set that earns the version bump. Everything else is post-3.0.

| Area | Items | Why in 3.0 |
|------|-------|-----------|
| **Object model** | SE return-type migration ([ADR-0002](adr/0002-return-type-migration.md)); return-type consistency (#71); `getGEO` extension policy ([ADR-0003](adr/0003-getgeo-extension-policy.md)) | The headline. (The premature NEWS SummarizedExperiment claim has been removed; `getGEO` still returns `ExpressionSet` until this lands.) |
| **Single-cell core** | SC architecture ADR (0004); SC1 manifest → SC2 10x → SC3 formats → SC4 combine | The other headline; the most-requested modern data class. |
| **Correctness floor** | #58/#154 `findFirstEntity` hardening; accession validation; #60 parseCharacteristics; #21 GDS NA ids; #131 URL join; #147 timeout floor | A 3.0 with known crash-on-input bugs is not paper-ready. |
| **Foundations** | offline test fixtures (XL); structured `rlang` conditions; BiocFileCache caching | Prerequisite for stable CI, coverage, and the paper's reproducibility story. |
| **Docs** | fill class/accessor roxygen (#103); restructure main vignette + complete SC vignette (#156); rewrite DESCRIPTION/biocViews; correct NEWS | **Hard prerequisite for the paper** — cannot cite a tool whose class docs are empty. |

### Post-3.0 (after the milestone, before or alongside the paper)

Triage automation (Layers 1–3), CI expansion (coverage/BiocCheck/lint/integration), retry/resume + parallel download, metadata-only fetch, tidy output, SOFT `series_table` (#80), token auth (#154), SC5–SC9 (lazy assays, spatial), remaining bug backlog (#98, #148, #14).

### Paper (capstone — strictly downstream of the 3.0 batch + docs)

A citable update is warranted: ~18 years since the original (Davis & Meltzer, *Bioinformatics* 2007), and GEO itself transformed (RNA-seq, single-cell, spatial). **Write about shipped reality, not this roadmap.**

- **Venue:** bioRxiv preprint first (priority + README/CITATION link) → **Bioinformatics Application Note** or **F1000Research** "GEOquery 3.0:" update. This sidesteps JOSS's prior-publication concern (the 2007 note already exists); use JOSS only to emphasize the software-engineering modernization.
- **Material already in hand:** the ADRs + this ROADMAP + NEWS history are the "design decisions / what changed since 2007" section. Keep writing ADRs through the 3.0 batch — they are the methods notes.
- **Needs before submission:** the 3.0 batch shipped, docs fixed (#103/#156), one or two worked examples / light benchmarks on real GSEs, updated `CITATION.cff`.

### Sequence

1. Land the 3.0 batch (object model + SC core + correctness floor + foundations + docs).
2. Cut **GEOquery 3.0.0** at the Bioconductor release boundary.
3. Worked examples / validation on real GSEs.
4. bioRxiv preprint → App Note; update `CITATION.cff`.

> **Failure mode to avoid:** announcing 3.0 + paper on aspiration, then spending a year catching the code up to the abstract. Code and docs ship first; the paper describes what exists.

---

## Legend

**Priority** — P0: critical / blocking or actively breaking users · P1: high leverage, do soon · P2: valuable, schedule deliberately · P3: nice-to-have / opportunistic.

**Effort** — S: < ~half day, single function · M: ~1–3 days · L: ~1 week · XL: multi-week / cross-cutting.

**Tags** — `bug` correctness · `feature` new capability · `ux` ergonomics/API · `tech-debt` internal cleanup · `docs` · `ci-cd` · `automation`.

Issue links use GitHub `#NNN` references. Items with no issue are maintainer-identified.

---

## Issue staleness & feasibility triage

The 17 open issues span 2015–2026. Several predate the 2024 httr2 migration and the GSEMatrix/SummarizedExperiment changes, so they may be fixed, mis-framed, or no longer reproducible. **Verify before scheduling work.** This table is the gate that feeds everything below — do not build on an unconfirmed issue.

**Verdicts:** `valid` confirmed against current code · `needs-repro` plausibly fixed/changed since filing, re-run before acting · `re-scope` real underlying bug but the issue as written is wrong/conflated · `close` recommend closing (not a GEOquery bug, obsolete, or wontfix).

| Issue | Age | Verdict | Reasoning |
|-------|-----|---------|-----------|
| #158 single-cell 10x | 2026-01 | **valid** | Fresh, clear scope. Real gap (no mtx/h5/h5ad reader). |
| #156 vignette overhaul | 2025-07 | **valid** | Docs still history-heavy, no quick-start. Confirmed by docs audit. |
| #154 private/token data | 2024-09 | **valid** (niche) | Token URLs hit public path → 404. Real, but security-sensitive + not CI-reproducible. Low priority. |
| #148 encoding + GSE233761 | 2024-03 | **needs-repro / re-scope** | Two problems conflated (unused `encoding` arg + matrix parse error). Predates httr2; re-test GSE233761 against current code, then split. |
| #147 timeout ignored | 2023-11 | **re-scope** | Original symptom partly addressed by httr2 migration, BUT the `max(getOption("timeout"),120)` floor (getGEOfile.R:136) is a real live bug. Keep, narrowed to the floor. |
| #133 `app$vspace` error | 2023-01 | **close** | `cli`/pkg-env error ("attempt to apply non-function"), not GEOquery code. Stale dependency in reporter's env. Not reproducible. Close with needs-info/wontfix. |
| #131 `fs::path` | 2022-06 | **re-scope** | Underlying `//`-in-URL bug is real, but the proposed fix is **wrong** — `fs::path` mangles `https://` → `https:/`. Close the original framing; track the real fix (Quick win 4, `url_join`/`req_url_path_append`) as a fresh scoped issue. |
| #103 docs missing | 2020-09 | **valid** | Class/accessor roxygen still near-empty. Confirmed. |
| #98 varMetadata (GSE136400) | 2020-01 | **needs-repro** | Predates current parseGSEMatrix. Re-run GSE136400; the failing construction path (parseGEO.R:603) still exists, so likely valid — confirm. |
| #80 series_table blocks | 2019-01 | **valid** | `!series_table_begin/_end` still unhandled. Confirmed. |
| #71 DESCRIPTION outdated | 2018-06 | **valid** (trivial) | DESCRIPTION still microarray-only; GSEMatrix is now default. Cheap doc fix. |
| #68 quiet arg | 2018-03 | **valid** (trivial) | No quiet control still. Additive. |
| #60 parseCharacteristics=FALSE | 2017-12 | **valid** | Forwarding bug **confirmed** in code (param dropped at parseGEO.R:452). Real. |
| #58 non-public infinite loop | 2017-12 | **valid** (re-verify accession) | `findFirstEntity` `while(TRUE)` + length>1 `grepl` hazard still present. Hardening valid regardless of whether GSM1062236 still behaves identically. |
| #21 NA in GDS ids | 2015-09 | **valid** (trivial) | GDS2MA.R:43 still unguarded; GSEMatrix path already guards. Confirmed. |
| #18 getGPL=FALSE local | 2015-09 | **needs-repro** | Code analysis suggests the flag *does* forward now; symptom may be `findFirstEntity` misclassification or already fixed. Re-run before acting. |
| #14 `amount="quick"` | 2015-09 | **needs-repro / likely close** | Depends on the `acc.cgi` "quick" SOFT view still existing/behaving as in 2015. Verify GEO still serves it; if obsolete, wontfix. Lowest priority. |

**Recommended closes now:** #133 (not our bug), #131 (wrong fix — re-file the real one). **Re-verify before any work:** #148, #147, #98, #18, #14. **Confirmed-real, safe to schedule:** #60, #58, #21, #71, #68, #103, #80, #156, #158, #154.

> This triage feeds the automation rubric: only **confirmed-valid + reproducible** issues should ever reach the issue→PR path. `needs-repro` issues go through Layer 2 (repro harness) first; `close`/`re-scope` issues are never auto-dispatched.

---

## Quick wins (P0/P1, effort S)

The highest-leverage small items. Each is single-function or single-file, has a concrete reproducer, and unblocks user-facing pain.

| # | Item | Pri | Eff | Issues | Problem → Approach |
|---|------|-----|-----|--------|--------------------|
| 1 | **Thread `parseCharacteristics` end-to-end** | P1 | S | #60 | `getGEO(parseCharacteristics=FALSE)` is a no-op: `getAndParseGSEMatrices` (R/parseGEO.R:432) accepts the flag but the `parseGSEMatrix` call (R/parseGEO.R:452–453) drops it, and the local-file path (R/getGEO.R:137 → `parseGEO`) never threads it. GSE23397 still crashes `tidyr::spread`. → Add `parseCharacteristics = parseCharacteristics` to the `parseGSEMatrix` call; add the param to `parseGEO()`'s signature and forward in the `0` branch; pass it from `getGEO()` at line 137. Add a local-fixture regression test asserting no characteristics columns are split. |
| 2 | **Harden `findFirstEntity` against private/HTML pages & length>1 grepl** | P0/P1 | S | #58, #154 | Private accession (GSM1062236) returns an HTML error page; `findFirstEntity` (R/parseGEO.R:217–242) loops `while(TRUE)` and the `if (!grepl(" = ", checkline))` test errors when `checkline` length>1 (R≥4.2). Result: unkillable hang or cryptic mis-parse. → Make the test scalar (`length(checkline)>0 && !any(grepl(" = ", checkline))`); add an early HTML/error sniff (`<html`/`<!DOCTYPE`/`not currently public`) that `abort()`s with a clear "accession may be private/embargoed/not public" message; bound the loop with a max-lines safeguard. |
| 3 | **Sanitize `NA` row names in `GDS2eSet`** | P2→P1 | S | #21 | `rownames(expr) <- as.character(Table(GDS)$ID_REF)` (R/GDS2MA.R:43) throws on NA ID_REF (GDS3666); the GSEMatrix path already guards this (parseGEO.R:596). → `idref[is.na(idref) \| idref==""] <- "NA"; make.unique(idref)`, apply same to the featureData match (line 55), add an NA-ID_REF fixture test. (Promoted because it is trivial and fully blocks a documented accession.) |
| 4 | **Fix `file.path`-in-URL double-slash; add `url_join` helper** | P1 | S | #131 | `file.path(url, fnames)` (getGEOSuppFiles.R:142) and `file.path(url,"filelist.txt")` (:166) mangle the `https://` scheme. The proposed `fs::path` fix in #131 is **wrong** (collapses `https://`→`https:/`). → Prefer `httr2::req_url_path_append()` (already used at :133); add a small `url_join()` for plain-string builds. Non-network unit test asserting no `//` after scheme. |
| 5 | **`quiet` argument for `getGEOSuppFiles`** | P1/P3 | S | #68 | No way to silence "Using locally cached version…" / "No supplemental files found"; `downloadFile` already hardcodes `quiet=TRUE`. → Add `quiet = getOption("GEOquery.quiet", FALSE)`, wrap `message()` calls and httr2 progress in `if(!quiet)`. Mechanical, no behavior change for existing callers. |
| 6 | **Remove the 120s timeout floor** | P1/P2 | S | #147 | `timeout_seconds <- max(getOption("timeout"),120)` (R/getGEOfile.R:136) silently overrides any user-set lower timeout and overloads base R's `timeout`. → Add a dedicated `GEOquery.download.timeout` option (default 300) + a `timeout` arg with **no** `max()` floor; replace the hardcoded `req_timeout(15)` in rnaseq.R for consistency. |
| 7 | **Add coverage workflow + Codecov badge** | P1 | S | — | `covr` sits unused in Suggests; README shows only R-CMD-check badges. → Add `.github/workflows/test-coverage.yaml` (r-lib recipe, `codecov/codecov-action@v5` with token) and the badge. **Sequence after fixtures** (item below) or run in the integration job — covr replays the same live-network tests and will be near-zero/flaky otherwise. |

> Items 1, 2, 4, 5 are also the first batch for the issue→PR automation pipeline (see Triage section). Treat them as both manual quick wins and pipeline shakedown candidates.

---

## Bugs & Correctness

| Title | Pri | Eff | Issues | Problem → Approach |
|-------|-----|-----|--------|--------------------|
| `getGPL=FALSE` ignored for local series_matrix | P2 | S | #18 | On the local path, `getGPL` does forward into the `0` branch (parseGEO.R:46) and the only GPL trigger is `if(getGPL)` at parseGSEMatrix:584 — so the symptom most likely stems from `findFirstEntity` misclassifying the local file (not taking the `0` branch). → Add a regression test asserting no network GPL fetch and 0 featureData columns; if detection misfires, fix `findFirstEntity` classification. Shares root with item 1 (parseCharacteristics also unforwarded here). |
| `AnnotatedDataFrame` varMetadata mismatch (GSE136400) | P2 | M | #98 | `new("AnnotatedDataFrame", data=dat, varMetadata=vmd)` (parseGEO.R:603–605) fails when vmd row count ≠ `ncol(dat)` after `make.unique`, or vmd lacks `labelDescription`. → Reindex vmd to `colnames(dat)`, ensure a `labelDescription` column (mirror GDS2eSet:58–65), wrap construction in validation with a clearer error. |
| `encoding` parameter unused; GSE233761 errors | P2 | M | #148 | `getGEO()` has no `encoding` arg (R/getGEO.R:117–118); all reads use `data.table::fread` defaults, so non-UTF-8 matrices mis-parse and downstream `t()`/`spread` error. → Either thread `encoding` into the `fread` calls (parseGSEMatrix:481/491/494) **or** remove it and emit a `deprecate_warn`/typed warning when supplied — never silently drop a user arg. Wrap the table read in `tryCatch` re-aborting with a `geoquery_parse_error` that hints "retry with GSEMatrix=FALSE". Reproduce GSE233761 with a saved fixture. |
| `downloadFile` streams into memory; large GSE aborts | P2 | S | #147 | `resp_body_raw` buffers the whole response and a single `req_timeout` deadline kills large `family.soft.gz` mid-stream. → Stream to disk via `httr2::req_perform(path=destfile)` (as `getGEOSuppFiles` already does); document `req_timeout` as the overall deadline. Pairs with the retry/resume feature below. |
| `amount="quick"` SOFT parsing fails (invalid `n`) | P3 | S | #14 | Truncated `acc.cgi` SOFT view yields a negative `numberOfLines` in `fastTabRead` / `.parseGPLWithLimits` (readLines `n`). → Clamp `numberOfLines` to ≥0 (or −1 = read-all); guard the with-limits parsers against negatives; handle missing `table_end` gracefully. Add a saved truncated-SOFT fixture test. |
| `printHead` corrupted definition in classes.R | P2 | S | #103 | R/classes.R:22–39 repeats `printHead <- function(x) # …` ~18 times; it parses only because the trailing body wins. → Replace with one clean definition + the existing body. Pure readability/safety, no behavior change. |

---

## New Features

| Title | Pri | Eff | Issues | Problem → Approach |
|-------|-----|-----|--------|--------------------|
| **BiocFileCache-backed caching** | P1 | L | — | `tech-debt`. Caching is `file.exists()` against `tempdir()` (getGEOfile.R:118–124, getGEOSuppFiles.R:122–132): files vanish across sessions, no staleness/integrity check, partial downloads treated as cached, parallel sessions race. → Add `R/cache.R` with `geoCache()` returning a persistent `BiocFileCache` (default `tools::R_user_dir("GEOquery","cache")`); route `downloadFile`/`getGEOSuppFiles` through `bfcrpath`/`bfcadd` keyed on URL; add `clearGEOCache()`; gate behind an option for one cycle for back-compat. Add BiocFileCache to Imports. |
| **Retry/resume + configurable timeout download path** | P1 | M | #147, #68 | Combines the #147 fix with robustness. → Rework `downloadFile` to `req_perform(path=destfile)` + `req_retry(max_tries, is_transient)` + option-driven timeout (no floor); add optional Range-header resume via a `.part` file renamed on completion; apply the same wrapper in `getGEOSuppFiles`. |
| **Metadata-only fetch** | P1 | M | #80 | No fast path to sample/series metadata without pulling the full matrix/GPL. → Add `getGEOMetadata(GEO)` (or `metadataOnly=TRUE`) that streams header lines up to `!series_matrix_table_begin`, reusing the `!Sample_` extraction from commit 52fe93d; `getGPL` defaults FALSE in this mode. |
| **Single-cell & spatial support** | P1 | XL | #158, #80 | Most-requested modern data class; today an SC `getGEO()` returns a near-empty object. **Promoted to its own section** → see [Single-Cell & Spatial](#single-cell--spatial). |
| **Auto-decompression for `getGEOSuppFiles`** | P3 | S | #158, #68 | `.gz/.tar.gz` left compressed; the common `_RAW.tar` single-cell case needs extraction. → Add `decompress=FALSE`; when TRUE, `untar()` archives into per-archive subdirs and `R.utils::gunzip()` standalone `.gz`, returning extracted paths as new columns. R.utils already in Imports. |
| **Parse `!series_table_begin/_end` SOFT blocks** | P2 | M | #80 | Only `!series_matrix_table_begin/_end` is handled (parseGEO.R:485–486); SOFT-level series tables are silently dropped. → Extend the `parseGSE` SOFT path to read these blocks with `fread(..., fill=TRUE)` (mirroring the #162 handling), attach to the GSE object, document the accessor. |
| **Private/token authenticated access** | P2 | M | #154 | Reviewer `token=XXXX` URLs hit the public path → 404. → Add a `token` arg to `getGEO`/`getGEOfile`/`getGEOSuppFiles`, thread into httr2 as `req_url_query(token=...)`, honor `GEO_ACCESS_TOKEN` env fallback, add a token-aware `acc.cgi` URL branch. Document the reviewer workflow. **Security-sensitive — human-only path** (see Triage). |
| **Tidy/long data.frame output** | P2 | M | — | No first-class tidy path; S4 objects are unfamiliar to tidyverse users. → Add `getGEOTidy()`/`as_tibble` that pivots the assay with `tidyr::pivot_longer` and left-joins phenoData by sample id. tidyr/dplyr already in Imports. |
| **Parallel/async multi-file download** | P3 | M | — | `getGEOSuppFiles` and multi-platform matrix fetches loop sequentially. → Use `httr2::req_perform_parallel()` with a modest `max_concurrent` (default 4, NCBI courtesy) + progress bar; keep cache-skip before scheduling. |

---

## Single-Cell & Spatial

> **The reframe that makes this section necessary:** for single-cell data the series matrix is the *wrong abstraction*. GEOquery's core model (assay matrix + phenoData → ExpressionSet) assumes data lives in the series matrix. SC data almost never does — it lives in **supplementary files**, and the series matrix is typically empty or metadata-only. So `getGEO("GSExxxxx")` on an SC study returns a near-useless object today. The fix is a **supplementary-files-first SC entry point**, not patching the matrix path. Every item below builds on that premise. Target output classes: `SingleCellExperiment`, and `SpatialExperiment` for spatial.

### Dependency architecture — DECIDED ([ADR-0004](adr/0004-single-cell-architecture.md))

Readers: **TENxIO** for 10x Matrix-Market (`.mtx`) and CellRanger HDF5 (`.h5`), **anndataR** for AnnData (`.h5ad`) — both Bioconductor, used for their `SingleCellExperiment` assemblers. **In-package**, in `Suggests`, behind `requireNamespace()` guards. `DropletUtils` (for `emptyDrops`/advanced CellRanger) and `BPCells` (on-disk/lazy, SC5) stay optional. Rejected as defaults: `DropletUtils` (heavy compiled chain — it broke Windows CI, #166), `zellkonverter` (bundles Python via basilisk), `Seurat` (weight/ecosystem), and hand-rolled `Matrix::readMM` (prefer Bioc assemblers). Prefer Bioconductor readers wherever one exists.

### Roadmap

| # | Title | Pri | Eff | Issues | Problem → Approach |
|---|-------|-----|-----|--------|--------------------|
| SC1 | **SC manifest / inventory primitive** *(foundational)* | P1 | M | #158 | No way to see what an SC GSE contains before a multi-GB download. → `geoSingleCellManifest(GSE)`: list supp files via `getGEOSeriesFileListing()` (no download), detect format + role per file, **group 10x MEX triplets by GSM prefix**, peek inside `_RAW.tar` listings; return a tibble `{gsm, format, role, url, size_bytes}`. Every item below consumes this. |
| SC2 | **10x MEX triplet assembly** | P1 | M | #158 | The real pain isn't reading — it's that GEO ships `GSMxxx_matrix.mtx.gz`/`_barcodes.tsv.gz`/`_features.tsv.gz` flat or tarred, while a reader wants a *directory per sample with canonical names*. → Auto-group + symlink/rename into the expected layout, then import via **TENxIO** ([ADR-0004](adr/0004-single-cell-architecture.md)) → `SingleCellExperiment`. ~80% of "handle 10x". |
| SC3 | **Multi-format reader dispatch** | P1 | L | #158 | Dispatch keyed on detected format → `SingleCellExperiment`, all `requireNamespace()`-guarded ([ADR-0004](adr/0004-single-cell-architecture.md)): `.mtx` + CellRanger `.h5` → **TENxIO**; `.h5ad`/AnnData → **anndataR** (native R, no Python). `DropletUtils` optional for `emptyDrops`/advanced CellRanger. h5ad is very common on recent GEO and currently forces hand-rolling. |
| SC4 | **Combined SCE + sample→cell metadata broadcast** | P1 | M | #158 | Per-GSM SCEs need merging, and cell barcodes don't map to GSM directly. → Merge per-sample SCEs into one with a `gsm`/`sample` colData column, **broadcasting GSM characteristics to every cell** (GEOquery already parses GSM metadata — unique advantage); handle feature-set union/intersection mismatch across runs. |
| SC5 | **On-disk / lazy assays for scale** | P1 | M | — | SC matrices blow past memory. → `delayed=TRUE` → HDF5Array/DelayedArray-backed assays, or **BPCells** (`open_matrix_10x_hdf5`, on-disk bit-packed, ~70× less memory) as an opt-in backend. Pairs hard with the **BiocFileCache** item (huge downloads must persist across sessions). Without this, large GSEs OOM. |
| SC6 | **Cell-level metadata join** | P2 | M | #80 | Many GSEs ship a separate `*_metadata.txt.gz` / barcode→cluster annotation (the #80 `!series_table` is one form). → Detect and join to `colData` **by barcode**. Turns raw counts into an analysis-ready annotated object. Overlaps the SOFT `series_table` parser in New Features. |
| SC7 | **Auto-detect + nudge in `getGEO()`** | P2 | S | #158 | SC `getGEO()` silently returns an empty object → user confusion. → Heuristic (mtx/h5ad/loom supp files present, or `scRNA`/library-strategy in metadata) emits `message("looks single-cell; use getGEOSingleCell()")`. Cheap, high UX payoff. |
| SC8 | **Pre-download size guard / selective fetch** | P2 | S | #158 | SC supp files reach tens of GB. → Surface total size from the SC1 manifest; let users filter to specific GSMs before download. Trivial once SC1 exists. |
| SC9 | **Spatial transcriptomics (Visium)** *(forward-looking)* | P3 | L | — | Spatial is growing fast on GEO and has no path today. → Same supp-files-first pattern, output `SpatialExperiment`: detect spaceranger output (`tissue_positions`, `scalefactors`) and read via `read10xVisium`. |

**Sequencing within SC:** ADR → SC1 (manifest) → SC2 (10x) + SC7 (nudge) → SC3 (formats) + SC4 (combine) → SC5 (lazy, with BiocFileCache) → SC6/SC8 → SC9 (spatial). SC1 is the linchpin; SC5 depends on the BiocFileCache caching item.

---

## UX & API

| Title | Pri | Eff | Issues | Problem → Approach |
|-------|-----|-----|--------|--------------------|
| **Accession-type validation with actionable errors** | P0 | M | #58, #148 | `geotype <- toupper(substr(GEO,1,3))` (getGEO.R:130, getGEOfile.R:46) is never validated; unknown prefixes leave `destfile`/`myurl` undefined ("object destfile not found") and feed the #58 hang. → Add `R/validate.R::validate_geo_accession()` with regex `^(GDS\|GSE\|GSM\|GPL)\d+$`, `rlang::abort(class="geoquery_bad_accession")`; call at the top of `getGEO`/`getGEOfile`. Combine with the `findFirstEntity` hardening (Bugs item 2). |
| **Structured condition classes (`rlang::abort`)** | P1 | L | #133, #58 | Every error is a bare `stop()`; `downloadFile` does `message(e)` then a generic `stop("Failed to download…")` (getGEOfile.R:153–167), discarding the httr2 status/curl condition (#133). → Add `rlang` to Imports; define `geoquery_error` parent + `geoquery_download_error` (carry status/url), `geoquery_parse_error` (fname), `geoquery_bad_accession`, `geoquery_private_accession`. Preserve the chain via `parent=e`. Document classes for downstream `tryCatch`. Foundational for items 2, #148, #133. |
| **Consistent `quiet`/verbose control** | P1 | M | #68 | Bare `message()` everywhere; `downloadFile`'s `quiet` not threaded up; `getAndParseGSEMatrices` unconditionally messages (parseGEO.R:439–443). → Single `quiet = getOption("GEOquery.quiet", FALSE)` on all download/parse fns + internal `inform(quiet, …)` wrapper; register the option in `R/zzz.R`. Supersedes the standalone `getGEOSuppFiles` quiet quick win. |
| **Predictable, documented return type** | P1 | L | #71 | `getGEO` returns a named list for GSE+matrix but bare S4 otherwise; callers must defensively `is.list()` + `[[1]]`. → **Always** return a list for GSE; add `simplify=FALSE` for single-element convenience; rewrite `@return` and DESCRIPTION to lead with the GSEMatrix default; add a migration `@section`. |
| **Deliver the SummarizedExperiment migration** | P1 | L | #71, #168 | `tech-debt`. The series-matrix path still returns `ExpressionSet` (`parseGSEMatrix`, parseGEO.R:620); the premature NEWS SE claim was removed in the version reconciliation. → Build SE directly (or `eset_to_se()`), gate behind `getGEO(..., returnType = c("ExpressionSet","SummarizedExperiment"))` defaulting to ExpressionSet for one release with a deprecation warning, then flip. Export `as_SummarizedExperiment()`. |
| **Validate `destdir`/`baseDir`/flag inputs** | P2 | M | #58, #68 | Paths/flags used directly in `file.path`/`dir.create`; `suppressWarnings(dir.create())` (getGEOSuppFiles.R:111) hides un-creatable dirs. → `R/validate.R::check_args()` with `geoquery_bad_input`: assert length-1 char paths (create + check return), logical scalars for flags, and a length-2 increasing positive-integer `GSElimits`. |

---

## Documentation

| Title | Pri | Eff | Issues | Problem → Approach |
|-------|-----|-----|--------|--------------------|
| Fill empty class roxygen (6 classes) | P1 | M | #103 | Every class block in R/GEOquery-package.R is one-line boilerplate with misleading `new(...)` advice. → Document each slot (from R/classes.R), replace `new(...)` text with "returned by `getGEO()` when `GSEMatrix=FALSE`", add `@examples` (`Meta()`/`Table()`), cross-link accessors. |
| Document GEOData accessor generics | P1 | M | #103 | Accessor block (lines 3–12) only points elsewhere; no `@param`/`@return`/`@examples` for `Meta`/`Table`/`Columns`/etc. → Add per-generic descriptions, params, returns, and a runnable example; regenerate `man/`. |
| Rewrite DESCRIPTION + biocViews | P1 | S | #71 | Description (line 35) still says microarray-only "bridge"; biocViews omit RNASeq/SingleCell. → Multi-sentence description covering Series Matrix→ExpressionSet, SOFT→S4, RNA-seq counts, search, supp/single-cell; add `GeneExpression, Transcriptomics, RNASeq, Sequencing, SingleCell, ThirdPartyClient`. |
| Correct `getGEO` return docs | P1 | S | #71, #156 | Roxygen leads with SOFT as primary though `GSEMatrix=TRUE` is default. → Lead Details/@return with the default behavior (list of `ExpressionSet`). (The NEWS version/SE mismatch was already fixed in the version reconciliation.) |
| Restructure main vignette `GEOquery.qmd` | P1 | L | #156 | History-heavy, no quick-start, never shows `exprs`/`pData`/`fData`. → Add a "Quick start" then task sections (download+access, GPL annotation, supp files, RNA-seq counts, search, GDS conversion, SOFT parsing); move history to a collapsible appendix. |
| Complete single-cell vignette | P1 | M | #158, #80, #156 | `single-cell.qmd` is a skeleton with an empty h5ad section and a hardcoded `/Users/davsean/...` path (line 100). → Add intro on GEO single-cell conventions, fill the h5ad section, replace the absolute path with `s$filepath`-derived value, add prose + expected output. Note manual workflow pending #158. |
| "Migrating from ExpressionSet" note | P2 | S | #156, #71 | Three object models with no map between them. → Add a concise table mapping task → ExpressionSet (`exprs`/`pData`/`fData`) vs SOFT (`Table`/`Meta`/`Columns`) vs SE (`assay`/`colData`/`rowData`). |
| Organize pkgdown reference + navbar | P2 | S | #156, #103 | `_pkgdown.yml` is 4 lines; reference page is a flat alphabetical dump. → Add grouped `reference:` sections + articles navbar; mark internal helpers `@keywords internal`. |
| Update README usage | P3 | S | #156 | Links to `GEOquery.Rmd` on `master`; vignette is now `.qmd`, default branch `devel`; no inline example. → Fix link to the pkgdown article, add a minimal `getGEO` example, update `master`→`devel`. |

---

## CI/CD & Engineering

| Title | Pri | Eff | Issues | Problem → Approach |
|-------|-----|-----|--------|--------------------|
| **Fixturize network tests (offline + deterministic)** | P1 | XL | — | `tech-debt`. All 8 testthat files hit **live** NCBI with zero `skip_on_cran`/`skip_if_offline` and no fixtures → slow, flaky, non-deterministic; covr unreliable. → Add `httptest2` to Suggests; record fixtures once with `capture_requests({...})` for the tiny accessions already used (GSE11413/11595/34145); wrap test bodies in `with_mock_dir(...)`. Gate FTP/RNA-seq scraping paths (bypass httr2) with `skip_if_offline()`/`skip_on_cran()`; keep heavy assertions in the integration job. **This is the linchpin** — it unblocks reliable coverage, lint, BiocCheck, and safe PR automation. |
| Coverage workflow + badge | P1 | S | — | (See Quick win 7.) Sequence after fixtures. |
| BiocCheck workflow | P1 | M | — | Bioc package but CI runs only `R CMD check`; #71/#103 surface only on the Bioc build report. → Add `.github/workflows/bioc-check.yaml`, ideally container-based (`bioconductor/bioconductor_docker:devel`) running `BiocCheck::BiocCheck()`. Downgrade WARNINGs initially given known doc issues, then tighten to fail-on-warning. |
| lintr workflow | P2 | S | — | `.lintr` committed but unenforced. → Add `.github/workflows/lint.yaml` (r-lib template) with inline PR annotations; `continue-on-error: true` for the first PRs to surface backlog, then make required. |
| Scheduled integration workflow (live NCBI) | P2 | M | #21, #58, #18 | Offline tests can't catch NCBI format drift (the class behind #162/#148/#60). → `.github/workflows/integration.yaml` on weekly cron + `workflow_dispatch`, `GEOQUERY_INTEGRATION=true` env to bypass mocks; on failure open/append a tracking issue. Keeps PR CI fast while catching upstream changes. |
| Concurrency + trigger hygiene | P2 | S | — | `pkgdown.yaml` rebuilds on every branch push; no `concurrency` groups → redundant overlapping runs. → Add `concurrency: {group, cancel-in-progress: true}` to all workflows; scope pkgdown to `[main, devel]` + doc paths; add `use-public-rspm: true` to the style job; bump `checkout@v3`→`@v4` in pr-commands.yaml. |
| Dependabot for Actions | P3 | S | — | No `dependabot.yml`; action pins already drift (`checkout@v3` vs `@v4`). → Add `.github/dependabot.yml` (`github-actions`, weekly, grouped minor/patch). |
| Pre-commit hooks | P3 | M | #103, #71 | Recurring doc/NAMESPACE drift; the `/document` and `/style` commands exist precisely because of it. → Add `.pre-commit-config.yaml` (lorenzwalthert hooks: roxygenize, styler, lintr, use-tidy-description) + a CI mirror; document setup in CONTRIBUTING. Shifts the manual slash-commands left to commit time. |

---

## Issue Triage & PR Automation  *(maintainer focus area — most detailed)*

Today there is exactly one automation: `pr-commands.yaml` (`/document`, `/style`, OWNER/MEMBER-gated on PRs). The 17 open issues land unlabeled and untriaged; the single maintainer reads each one and hand-runs `getGEO` to confirm accession-specific bugs. The design below is a **three-layer system**, each layer independently useful, with security and human-gating baked in. The maintainer already uses agentic coding agents (PR #162 was a `copilot/fix-getgeo-error` branch), so the agent path is realistic — but it must be gated.

### Layered design

**Layer 1 — Auto-triage on issue open** (`P1`, `M`, `.github/workflows/issue-triage.yaml`)
Trigger `issues: [opened, edited, reopened]`, `permissions: {issues: write, contents: read}`. One job, two stages:
- **Deterministic pre-pass** (bash/Rscript, **source of truth**): grep the body for accession regex `\b(GSE|GSM|GPL|GDS)[0-9]{1,9}\b` (case-insensitive, dedup, cap ~10); keyword→label heuristic (`error|traceback|crash|infinite loop`→`bug`; `vignette|docs|documentation`→`docs`; `how do I|support`→`question`; else `enhancement`). Apply the `has-accession` label and list detected accessions in the comment to feed Layer 2.
- **Claude pass** (advisory) via `anthropics/claude-code-action@v1` in **text-only mode** (no shell tools), `ANTHROPIC_API_KEY` secret: produce labels from a **fixed allowlist** `{bug, enhancement, documentation, question, needs-repro, has-accession}` plus a 3–5 line triage comment. Parse its JSON, **intersect with the allowlist** before applying, so a prompt-injected body cannot create arbitrary labels.
- **Security:** the issue body is attacker-controlled — pass it via env var / stdin file, **never** interpolate into a shell command.

**Layer 2 — Reproduction harness** (`P1`, `L`, `.github/workflows/repro-harness.yaml`)
Trigger `issues: [labeled]` gated to `needs-repro`, plus `workflow_dispatch` with `accession`+`issue_number` inputs. Runs in `container: bioconductor/bioconductor_docker:devel`. Reads the accession **only** from the deterministic extractor / dispatch input (never re-parses arbitrary body text), installs the package, runs `tryCatch(withCallingHandlers(getGEO(acc, GSEMatrix=TRUE), warning=...), error=…)` capturing stdout/stderr/warnings/`sessionInfo()`/`packageVersion()`. Hard `timeout-minutes: 20` **and** an R-side `withTimeout`/`req_timeout` so the #58/#147 hangs can't pin a runner. Truncate to ~64KB, upload full log as artifact, post a collapsible PASS/FAIL comment. Re-run weekly over all open `has-accession` issues to detect "fixed upstream / now reproduces" drift. **Security:** downloads only NCBI accessions, runs no issue-supplied code, uses `GITHUB_TOKEN` only — safe for non-collaborator issues.

**Layer 3 — Issue→PR dispatch** (`P1`, `L`, `.github/workflows/issue-to-pr.yaml`)
Trigger **only** `issues: [labeled]` with `if: github.event.label.name == 'agent-ready'` — a **human-applied** label is the sole dispatch gate (Layer 1 may *suggest* it in a comment but never apply it). Enforce `author_association ∈ {OWNER, MEMBER}` (mirror pr-commands.yaml). `permissions: {contents: write, pull-requests: write, issues: write}`. Runs `claude-code-action@v1` in agent mode: branch `agent/issue-<n>-<slug>`, minimal change + targeted test, `devtools::document()` + `R CMD check`, open a **draft** PR ("Fixes #n", labels `agent-authored` + `needs-human-review`).
**Guardrails (workflow + CLAUDE.md):** (1) never auto-merge — draft PR, branch protection on `devel` requires the R-CMD-check matrix + human approval; (2) restrict editable paths to `R/ man/ tests/testthat/ NEWS DESCRIPTION` — **block `.github/workflows/**`** so the agent cannot rewrite its own guardrails or exfiltrate secrets; (3) require `has-accession` + a Layer-2 repro for any code-path fix; (4) run the full R-CMD-check matrix on the agent PR. **Security:** never expose `ANTHROPIC_API_KEY` on `pull_request` from forks — the agent runs in the trusted base-repo `issues` context only.

### Issue → PR decision rubric

A lightweight eligibility check comments *eligible / not-eligible* rather than auto-dispatching.

| | GOOD candidate | BAD candidate |
|---|----------------|---------------|
| Scope | Single function, < ~40 LOC delta | Multi-file; touches parser data model |
| Type safety | No S4 signature change; no change to `GSEMatrix` default or SE return contract | Changes class/return contract; new dependency |
| Testability | Behavior is one testable assertion | Ambiguous / under-specified |
| Evidence | Concrete failing accession or non-network repro exists | Docs-overhaul or design question |

### Current open issues as automation candidates

| Issue | Verdict | Why | Route |
|-------|---------|-----|-------|
| **#68** quiet arg | **GOOD — first** | Purely additive param in `getGEOSuppFiles`; `downloadFile` already `quiet=TRUE`. **No network needed to test.** | `agent-ready` |
| **#131** `//` in URL | **GOOD — first** | Bug literally at getGEOSuppFiles.R:142/:166; tiny `url_path_append` helper; **non-network** unit test. | `agent-ready` |
| **#60** parseCharacteristics | **GOOD** | One forwarding chain + a guard; concrete failing accession (GSE23397) for a Layer-2 repro. | `agent-ready` (after repro) |
| **#18** getGPL=FALSE local | **GOOD** | Single conditional short-circuit; targeted test against the accession. | `agent-ready` (after repro) |
| **#21** GDS NA ids | GOOD-ish | Trivial sanitize, but needs a live/fixture GDS3666 repro. | `agent-ready` (after fixture) |
| **#148** encoding + GSE233761 | **BAD** | Two intertwined problems + a parser data-model question — **split first**. | `needs-design` |
| **#133** `app$vspace` cli error | **BAD** | Environment/dependency-version bug, not deterministically reproducible. | human + Layer-2 repro |
| **#154** token auth | **BAD** | Auth + httr2 layer; needs a private accession + secret — **not CI-reproducible**, security-sensitive. | human only |
| **#156** vignette overhaul | **BAD** | Open-ended docs/UX, no acceptance criteria. | `needs-design` (split into scoped sub-issues) |
| **#158** single-cell 10x | **BAD** | New format parsers, new deps, multi-file, no single assertion. | `needs-design` / `help-wanted` |

> The eligibility check should auto-apply `needs-design`/`help-wanted` to BAD candidates with the reason commented, so the agent path is never pointed at them.

---

## Suggested sequencing

### Now (Phase 1 — foundations + bleeding fixes)
- **Bugs:** #58/#154 `findFirstEntity` hardening (P0) + accession validation (P0); #60 parseCharacteristics forwarding; #21 GDS NA ids; #131 URL join; #147 timeout floor.
- **Engineering:** **fixturize network tests (XL)** — the linchpin; structured `rlang` condition classes (foundational for many error-handling items).
- **Triage:** ship **Layer 1** (auto-triage) and **Layer 2** (repro harness) — both are read-only/safe and immediately reduce maintainer load.
- **Docs (cheap, high-trust):** correct `getGEO` return docs + NEWS mismatch; rewrite DESCRIPTION/biocViews.

### Next (Phase 2 — leverage on the foundation)
- **Engineering:** coverage + BiocCheck + lint workflows (now meaningful post-fixtures); scheduled integration job; concurrency/trigger hygiene.
- **Triage:** ship **Layer 3** (issue→PR) and run the **first batch (#68, #131, then #60, #18)** through it as a controlled shakedown.
- **Features/UX:** retry/resume download path (folds in #147); consistent `quiet`; predictable return type + the real SummarizedExperiment coercion (#71); BiocFileCache caching.
- **Docs:** restructure `GEOquery.qmd`; fill class/accessor roxygen (#103); pkgdown reference grouping.

### Later (Phase 3 — bigger bets)
- **Single-cell & spatial:** start with the dependency-architecture ADR, then **SC1 manifest → SC2 10x → SC3 formats / SC4 combine** (see [Single-Cell & Spatial](#single-cell--spatial)); SC7 nudge + auto-decompression are cheap and can land earlier. SC5 (lazy assays) waits on BiocFileCache; SC9 (spatial) is exploratory.
- **Features:** metadata-only fetch; token auth (#154, human-only); SOFT `series_table` blocks (#80); tidy output; parallel download.
- **Docs:** complete single-cell vignette; migration note; README example.
- **Engineering:** Dependabot; pre-commit hooks.
- **Bugs:** #98 varMetadata; #148 encoding (post-split); #14 `amount="quick"`.

### Dependency notes
- **Offline test fixtures must precede** aggressive PR automation (Layer 3) and the coverage/lint/BiocCheck signals — covr/lint replay the same live-network tests and will be flaky/near-zero otherwise; the agent path relies on deterministic tests to validate its draft PRs.
- **Layer 1's deterministic accession extractor feeds Layer 2**, and **Layer 2 repro results gate Layer 3** code-path fixes — build them in order.
- **Structured condition classes** should land before/with the error-surfacing items (#148, #133, accession validation) so all new errors are catchable from day one.
- The **SummarizedExperiment coercion (#71) and the doc/NEWS corrections must move together** to avoid re-introducing the code-vs-docs mismatch.
