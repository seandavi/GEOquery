# ADR-0005: Flip the getGEO default to SummarizedExperiment immediately

- Status: accepted
- Date: 2026-06-13
- Deciders: Sean Davis
- Amends: [[0002-return-type-migration]] (the migration *schedule* only)

## Context

[[0002-return-type-migration]] decided to migrate `getGEO()`'s GSE Series Matrix
return type from `ExpressionSet` to `SummarizedExperiment` through a **staged**
deprecation: ship the `returnType` argument defaulting to `ExpressionSet` for
one release cycle (with a "default will change" notice), then flip the default
in a later cycle.

Stage 1 shipped (the `returnType` argument, the `as_SummarizedExperiment()`
coercer, and the one-time notice). The maintainer has decided not to wait a
further release cycle for the flip — the modern object model is the headline of
this development series, and a long staged window mostly prolongs the
inconsistency between the docs/intent and the actual default.

## Decision

Flip the default **now**: `getGEO(returnType = c("SummarizedExperiment",
"ExpressionSet"))`, so the default is `SummarizedExperiment`.

- The `returnType` argument and `as_SummarizedExperiment()` from
  [[0002-return-type-migration]] are unchanged; only the default value of
  `returnType` changes.
- The one-time message is repurposed: instead of "the default will change", it
  now informs users that the default *has* changed and how to opt back into
  `ExpressionSet`.
- This is a genuine breaking change for code that uses `exprs()`/`pData()`/
  `fData()` on `getGEO()` results, and is recorded prominently under
  "Breaking changes" in NEWS.

## Consequences

- The package's default object model is now `SummarizedExperiment`, consistent
  with the single-cell (`SingleCellExperiment`) work in
  [[0004-single-cell-architecture]] and with current Bioconductor practice.
- Downstream code and reverse dependencies that assumed `ExpressionSet` break
  until updated; the escape hatch (`returnType = "ExpressionSet"`) and the
  coercer keep the migration cheap, and the change ships within a single devel
  series rather than spanning releases.
- Only the *schedule* of [[0002-return-type-migration]] is superseded; its
  technical decision (enum argument, coerce-at-boundary, `as_SummarizedExperiment()`)
  stands.

## Alternatives considered

- **Keep the staged schedule (flip next cycle):** safer for reverse deps, but
  the maintainer judged the prolonged inconsistency not worth it for a
  development series whose theme is the modern object model. Rejected.
