# ADR-0004: Single-cell reader architecture — TENxIO + anndataR, in-package, Suggests-guarded

- Status: accepted
- Date: 2026-06-13
- Deciders: Sean Davis

## Context

GEOquery will add single-cell support (roadmap SC items; entry point `getGEOSingleCell()`). Single-cell data on GEO lives in supplementary files, predominantly as 10x Matrix-Market triplets (`matrix.mtx` / `barcodes.tsv` / `features.tsv`), CellRanger HDF5 (`.h5`), and AnnData (`.h5ad`). We need readers that turn these into `SingleCellExperiment` objects.

The naive choice — `DropletUtils::read10xCounts()` — pulls a heavy compiled dependency chain (beachmat, BiocParallel, DelayedArray, edgeR, Rhdf5lib). That chain failed to build on the Windows CI runner and broke `R CMD check` (windows-only), traced and fixed in #166. We do not want the single-cell feature to drag that weight into the package's dependency surface, especially for the common mtx case.

A survey of the current ecosystem (June 2026):
- **TENxIO** (Bioconductor, waldronlab) — focused importer for 10x files; reads `.mtx` and `.h5`, returns `SingleCellExperiment` / `SummarizedExperiment` via Bioconductor's own assemblers. Much lighter than DropletUtils.
- **anndataR** (Bioconductor, scverse; v1.0.0, 2025) — native R reader for `.h5ad`/Zarr with **no Python**, converts to `SingleCellExperiment`. Supersedes zellkonverter (which bundles a Python/conda env via basilisk) for plain reading.
- **DropletUtils** — still best for `emptyDrops()` empty-droplet calling and CellRanger v2/v3 auto-detection; heavy.
- **BPCells** — fastest, on-disk, ~70× less memory; for very large data.
- Rejected for defaults: zellkonverter (Python weight), Seurat (large non-Bioconductor dependency tree, returns Seurat objects not SCE), and hand-rolling SCE from `Matrix::readMM` (we prefer Bioconductor's tested assemblers over maintaining our own).

This also resolves the in-package-vs-companion-package question deferred from [[0003-getgeo-extension-policy]].

## Decision

1. **Readers:** use **TENxIO** for 10x Matrix-Market (`.mtx`) and CellRanger HDF5 (`.h5`), and **anndataR** for AnnData (`.h5ad`). Both are Bioconductor packages; we rely on their assemblers to build `SingleCellExperiment` objects rather than constructing them ourselves.
2. **Placement:** single-cell support lives **in-package** (not a separate companion package). The reader packages go in **`Suggests`**, gated at call sites with `requireNamespace()` and an actionable error ("install TENxIO to read 10x files").
3. **Prefer Bioconductor:** when a Bioconductor reader exists for a format, use it over CRAN/Python alternatives, to inherit the ecosystem's assemblers, class conventions, and coercions.
4. **Optional enhancements (not default deps):** `DropletUtils` (for `emptyDrops`, advanced CellRanger handling) and `BPCells` (for on-disk/lazy large matrices, pairing with roadmap SC5) remain optional, also Suggests-guarded.
5. **Not used:** `zellkonverter` (Python/basilisk), `Seurat` (weight/ecosystem), and hand-rolled `Matrix::readMM` assembly.

## Consequences

- The common single-cell path stays light: no beachmat/DropletUtils compile chain in the default dependency surface, removing the cross-platform build fragility that caused the Windows CI failure (#166).
- Format coverage out of the box: `.mtx` + `.h5` (TENxIO) and `.h5ad` (anndataR) — the bulk of GEO single-cell submissions.
- Because readers are Suggests-guarded, the single-cell vignette and `getGEOSingleCell()` must `requireNamespace()` before use and degrade with a clear install message; the vignette must not unconditionally `library(DropletUtils)` (the #166 footgun).
- Advanced features (empty-droplet calling, out-of-memory analysis) are available but opt-in, keeping casual installs small.
- anndataR and TENxIO still depend on `rhdf5`/`Rhdf5lib` for HDF5; that is a lighter and more standard chain than DropletUtils', and only exercised when a user actually reads `.h5`/`.h5ad`.
- This is consistent with the extension policy in [[0003-getgeo-extension-policy]]: single-cell is a new verb (`getGEOSingleCell()`), not a `getGEO` flag.

## Alternatives considered

- **DropletUtils as the default reader:** most complete, but its heavy compiled chain broke Windows CI and bloats installs. Demoted to optional. Rejected as default.
- **Companion package `GEOquerySingleCell`:** keeps core leaner, but with Suggests-guarded readers the in-package weight is already negligible, and one install is friendlier to users. Rejected for now; revisit only if dependency or BiocCheck pressure grows.
- **zellkonverter for h5ad:** functional but bundles a Python/conda environment (basilisk) — exactly the install pain to avoid now that anndataR reads h5ad natively. Rejected.
- **Seurat readers / hand-rolled `Matrix::readMM`:** wrong ecosystem / reinventing assemblers we would have to maintain. Rejected in favour of Bioconductor importers.
