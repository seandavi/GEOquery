# ADR-0001: Record architecture decisions

- Status: accepted
- Date: 2026-06-13
- Deciders: Sean Davis

## Context

GEOquery has accumulated significant architectural decisions that are not obvious from the code alone — the split between the SOFT and Series Matrix parse paths, the `fastTabRead` speed/correctness tradeoff, regex-based dissection of `fread`-loaded character vectors, and the handling of malformed GEO records. These choices are currently captured only in scattered code comments and `NEWS.md` entries. New contributors (and future maintainers) repeatedly re-derive the reasoning.

## Decision

We will record significant architecture decisions as Architecture Decision Records (ADRs) in the `adr/` directory. Each ADR is a numbered Markdown file (`NNNN-title.md`) following `adr/template.md`. ADRs are immutable once accepted: to change a decision, write a new ADR that supersedes the old one and update the old one's status.

## Consequences

- Decision rationale is version-controlled alongside the code and reviewable in pull requests.
- The `adr/` directory is excluded from the R package build (`.Rbuildignore`) so it does not affect `R CMD check` or package size.
- Contributors must remember to write an ADR for non-trivial decisions; this is a process habit, not enforced by tooling.

## Alternatives considered

- **Wiki / external docs**: drifts from the code and is not reviewed with changes. Rejected.
- **Code comments only**: already in use and demonstrably insufficient for cross-cutting decisions. Rejected.
