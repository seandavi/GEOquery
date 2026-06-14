# ADR-0002: Migrate getGEO return type from ExpressionSet to SummarizedExperiment

- Status: accepted
- Date: 2026-06-13
- Deciders: Sean Davis

## Context

For the GSE Series Matrix path (the `GSEMatrix=TRUE` default), `getGEO()` constructs and returns Bioconductor `ExpressionSet` objects (`parseGSEMatrix()`, `new("ExpressionSet", ...)` at R/parseGEO.R:620). `ExpressionSet` is the legacy expression container; `SummarizedExperiment` is the modern Bioconductor standard and the required substrate for downstream classes (`SingleCellExperiment`, `SpatialExperiment`) we intend to support (see [[0003-getgeo-extension-policy]] and the single-cell roadmap).

Two complications:

1. **The NEWS file already claims the change shipped.** NEWS 2.99.0 (2024-10-01) states *"getGEO() now returns a list of SummarizedExperiment objects"*, but the code still returns `ExpressionSet`. The advertised breaking change was never actually implemented. This must be reconciled — users read NEWS.
2. **The two return types are not API-compatible.** `ExpressionSet` uses `exprs()`/`pData()`/`fData()`; `SummarizedExperiment` uses `assay()`/`colData()`/`rowData()`. A silent default flip breaks every downstream script and reverse-dependency. GEOquery sits high in the Bioconductor dependency graph, so a hard flip is hostile.

## Decision

We will migrate to `SummarizedExperiment` as the default return type **through a staged deprecation**, controlled by a single extensible argument rather than a boolean:

```r
getGEO(..., returnType = c("ExpressionSet", "SummarizedExperiment"))
```

- `match.arg()`, default = first element. The enum (not a `SummarizedExperiment = TRUE` boolean) is deliberate: it absorbs `"SingleCellExperiment"` later without a second conflicting flag.
- **Implementation: build once, coerce at the boundary.** Keep the existing `parseGSEMatrix()` assembly (phenoData/featureData/characteristics logic is correct) and coerce with `SummarizedExperiment::makeSummarizedExperimentFromExpressionSet()` when requested. One gated line, near-zero risk, no reimplementation. Also export a standalone `as_SummarizedExperiment()` coercer.

Migration schedule (one step per Bioconductor release cycle):

| Phase | Default | Opt-in | Signal |
|-------|---------|--------|--------|
| Now | `ExpressionSet` | `SummarizedExperiment` | deprecation message: default will change |
| Next | `SummarizedExperiment` | `ExpressionSet` | deprecation message on the ES opt-in |
| Later | `SummarizedExperiment` | `ExpressionSet` (legacy, retained) | — |

Immediately: **correct the NEWS 2.99.0 entry** to say `ExpressionSet` until the switch actually lands. The code and the docs move together (the doc/NEWS correction and the coercion must not drift apart again).

## Consequences

- No downstream breakage without a release of warning; reverse-deps get a cycle to adapt.
- A single argument serves both transition directions and future single-cell output, consistent with the extension policy in [[0003-getgeo-extension-policy]] ("output-shape args may join the signature; prefer a coercer").
- `returnType` is the one new flag we accept on `getGEO`'s signature; transport/IO concerns are pushed to options instead.
- Carrying both code paths (really: one path + a coercion) during the transition is a small maintenance cost.
- Per the policy ADR, this argument should still be considered against a coercer-first approach; we keep `as_SummarizedExperiment()` as the lower-friction escape hatch for users who do not want to thread the arg.

## Alternatives considered

- **Hard flip now (treat NEWS as the promise already made):** advertised ≠ shipped; reverse-deps never adapted to code that never returned SE. Rejected.
- **Boolean `SummarizedExperiment = TRUE`:** cannot grow to `SingleCellExperiment`; two competing booleans can conflict. Rejected in favour of the enum.
- **Native SE construction (drop ExpressionSet assembly):** rewrites correct, tested metadata logic for no functional gain during a transition that must keep ES available anyway. Rejected; coerce instead.
