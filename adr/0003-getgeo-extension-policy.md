# ADR-0003: getGEO extension policy — freeze the signature, route new capability deliberately

- Status: accepted
- Date: 2026-06-13
- Deciders: Sean Davis

## Context

`getGEO()` is the package's single high-level entry point with a large, long-standing user base. Its current signature is eight arguments: `GEO`, `filename`, `destdir`, `GSElimits`, `GSEMatrix`, `AnnotGPL`, `getGPL`, `parseCharacteristics`. It is tolerable today but shows three established overload smells:

1. **Mutually-exclusive entry modes.** `GEO` xor `filename` — fetch-remote vs read-local are two jobs behind one door.
2. **Flag-driven return type.** It returns `list<ExpressionSet>`, a single `ExpressionSet`, or an S4 `GSE`/`GSM`/`GPL`/`GDS` depending on flag values; callers defensively branch (issue #71).
3. **Arguments that silently no-op in some modes.** `parseCharacteristics`/`GSEMatrix`/`AnnotGPL` apply only to GSE; `GSElimits` only to SOFT-GSE. Issues #60 and #18 are bugs of exactly this shape — an argument ignored in a mode.

The improvement roadmap proposes substantial new capability (single-cell/spatial readers, metadata-only fetch, token auth, encoding, quiet, configurable timeout, decompression, `returnType`). Adding all of these as `getGEO` flags would push the function decisively into unmaintainable overload within one release. The overload is at the **signature** level, not the implementation level — internals are already factored (`getGEOfile` download, `parseGEO` dispatch, `parseGSEMatrix` assembly) — so the remedy is a gate on the signature, not a refactor.

## Decision

We adopt a standing policy governing what may be added to `getGEO()`:

1. **New data modalities and new verbs get their own functions, never a flag.** Single-cell → `getGEOSingleCell()`; metadata-only → `getGEOMetadata()`. No `singleCell=TRUE`, no `metadataOnly=TRUE` on `getGEO`.
2. **Cross-cutting IO/transport configuration goes to options (or env vars), not per-call flags.** `quiet` → `GEOquery.quiet`; download timeout → `GEOquery.download.timeout`; auth token → `GEO_ACCESS_TOKEN` / option. The `getGEO` signature stays frozen against transport concerns.
3. **Only genuine output-shape arguments may join the signature — and even then prefer a coercer.** `returnType` (see [[0002-return-type-migration]]) is the one accepted addition because it concerns output shape, not a new job; we still ship `as_SummarizedExperiment()` so users can avoid threading the arg.
4. **We do not split `getGEO`/`readGEO`** to resolve the `GEO`/`filename` smell. It is a real smell, but splitting is a breaking change to a heavily-used function and is not worth it. We document the constraint and live with it.

## Consequences

- `getGEO`'s signature stays stable; new functionality lands as discoverable sibling functions or options.
- Bugs of the "ignored-arg-in-mode" class (#60, #18) become less likely, because mode-specific behavior lives in mode-specific functions.
- Independently of arg count, we commit to fixing the return-type inconsistency (#71): always return a list for GSE plus a `simplify=` convenience. This reduces perceived overload more than trimming any single argument.
- A minor cost: capability is spread across several entry points, so documentation and the pkgdown reference must group them clearly (a roadmap docs item).
- The single-cell dependency-architecture decision (whether SC readers live in-package behind `Suggests` or in a companion package) is deferred to [[0004-single-cell-architecture]].

## Alternatives considered

- **Keep adding flags to `getGEO`:** simplest per-change, but the roadmap additions tip it into overload within a release. Rejected.
- **Aggressive refactor / split now:** breaks a large user base for an internal-cleanliness gain the factored internals do not require. Rejected.
