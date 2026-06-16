# ADR-0006: Seurat interoperability for single-cell data

- Status: accepted
- Date: 2026-06-13
- Deciders: Sean Davis
- Amends: [[0004-single-cell-architecture]] (the Seurat exclusion only)

## Context

[[0004-single-cell-architecture]] deliberately excluded Seurat: `Seurat .rds`
files were "not covered", and the readers returned only `SingleCellExperiment`.
The reasoning was to keep the dependency surface light and avoid the large
Seurat CRAN tree.

In practice two things make Seurat cheap to support:
1. The `SingleCellExperiment` <-> `Seurat` coercion is built into Seurat
   (`as.Seurat()` / `as.SingleCellExperiment()`) and is effectively lossless for
   our purposes.
2. GEO does ship single-cell data as saved Seurat objects in `.rds`
   supplementary files, which users currently cannot load through GEOquery.

So Seurat can be a *boundary coercion* rather than a new internal data model.

## Decision

Support Seurat at the edges, with `SingleCellExperiment` remaining the internal
lingua franca; `Seurat` is an optional (`Suggests`) dependency, `requireNamespace()`-guarded.

1. **Reading `.rds`** (`readGEOSingleCell`): a `.rds` is opaque from its name, so
   detect by class after `readRDS()` — a `SingleCellExperiment` is returned
   as-is, a `Seurat` object is coerced **to** `SingleCellExperiment` (via
   `Seurat::as.SingleCellExperiment`), and anything else is a clear error.
2. **Returning Seurat** (`readGEOSingleCell`/`getGEOSingleCell`): an
   `as = c("SingleCellExperiment", "Seurat")` argument coerces the
   `SingleCellExperiment` result **to** `Seurat` (via `Seurat::as.Seurat`) on
   output. Internally everything is read and combined as `SingleCellExperiment`;
   coercion happens only at the boundary.
3. **Documentation:** the single-cell article shows the coercion both ways for
   users who prefer to do it themselves.

## Consequences

- A saved Seurat object on GEO is now loadable; Seurat-first users can get a
  `Seurat` object directly; neither adds a hard dependency (Seurat stays
  `Suggests`, installed as a binary in the Bioconductor CI image).
- The read-`.rds`-Seurat → SCE and SCE → Seurat paths share Seurat's own
  coercions, so a Seurat round-trip is as lossless as Seurat itself makes it.
- loom and `_RAW.tar`-packaged / idiosyncratic layouts remain out of scope
  (unchanged from [[0004-single-cell-architecture]]).
- Only the Seurat exclusion of [[0004-single-cell-architecture]] is superseded;
  the TENxIO/anndataR reader choices and the in-package, Suggests-guarded
  placement stand.

## Alternatives considered

- **Keep Seurat fully out of scope (document coercion only):** lowest surface,
  but leaves saved-Seurat `.rds` files unreadable and forces every Seurat user
  to coerce manually. Rejected in favour of the thin boundary support.
- **Make Seurat an internal data model / hard dependency:** large, unnecessary;
  `SingleCellExperiment` already serves as the common representation. Rejected.
