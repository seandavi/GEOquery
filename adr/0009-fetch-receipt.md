# ADR-0009: Optional JSON receipts for GEO retrieval

- Status: proposed
- Date: 2026-09-27
- Deciders: Sean Davis (pending review); proposal by Cody Mitchell

## Context

Issue [#221](https://github.com/seandavi/GEOquery/issues/221) requests a small
retrieval manifest, including file hashes and available update dates. A GSE
can load multiple matrices and nested GPL files. Destdir reuse bypasses HTTP;
BiocFileCache may supply a file without its original response headers.
Related [#224](https://github.com/seandavi/GEOquery/issues/224) covers attaching
provenance to returned objects and remains separate.
The implementation is proposed in [#231](https://github.com/seandavi/GEOquery/pull/231).

ADR-0003 routes transport concerns through options rather than new getGEO
arguments. Existing checksum helpers and jsonlite suffice.

## Decision

An output path in `options(GEOquery.fetch.receipt = "new-receipt.json")` opts
getGEO and downloading getGEOSuppFiles calls into one JSON receipt per
successful outer call. NULL or FALSE disables the feature. The parent must
exist and be writable; an existing output path is refused. Listing-only
getGEOSuppFiles calls write no receipt. Returned objects remain unchanged.

A scoped collector includes nested downloads and local inputs. It is cleared
on exit, including errors. A nested getGEO call contributes to the outer
receipt rather than starting another file. A file reused during the same
call appears once, preserving its first download observation.

The top-level object has `schema_version: 1`, `operation` and a `resources`
array, including an empty array when no resources were retrieved. Each
resource has these fields, with JSON null for unknown values:

| Field | Meaning |
| --- | --- |
| accession | Parsed GEO accession when available, otherwise the retrieval URL's accession |
| requested_accession | Accession from the retrieval URL, or null for an unassociated local input |
| accession_mismatch | Whether parsed and requested accessions differ; null when either is unavailable |
| url | Sanitized request URL, or persistent cache key; null for local reuse |
| endpoint | That URL's origin and path, excluding its query |
| expected_url | For local reuse, the expected retrieval URL if known, without asserting file origin |
| file | Basename of the local input file |
| source | download, cache or local |
| bytes | Size of the bytes present on disk |
| md5 | MD5 of those exact bytes, using the existing checksum helper |
| http_last_modified | Last-Modified observed on this call's download, otherwise null |
| geo_last_update | Date from already parsed GEO metadata, otherwise null |
| geoquery_version | GEOquery version producing this receipt |
| timestamp | UTC time the stored bytes were observed, including on cache/local reuse |

The timestamp is not a claim about an earlier download. A cache key does not
establish current remote freshness. Checksums identify local bytes; they are
not comparison against a publisher-supplied checksum. Supplementary files
are not parsed, and no HTTP or GEO metadata requests are added.

Parsed accessions are compared after trimming surrounding whitespace and
normalizing case. A disagreement preserves both identities in the receipt
and raises a `geoquery_receipt_accession_mismatch` warning, inheriting from
`geoquery_warning`, once per resource. The condition carries `path`,
`requested_accession` and `observed_accession`. Returned objects are unchanged.
Normal R warning handling applies, including `options(warn = 2)`, which
aborts the call and prevents a completed receipt. With receipts disabled,
this identity check adds no warning or work.

Receipt URLs omit userinfo, fragments and all query parameters except GEO
routing fields `acc`, `targ`, `form`, `view` and `tool`. This excludes reviewer
tokens and arbitrary signed-URL credentials. Accessions and basenames may
still be sensitive for private records. Receipt sharing is a user decision.

A completed document is written to a temporary sibling file and renamed to
the requested path only after retrieval, parsing and return conversion
succeed. Receipt errors propagate. Existing paths, including collisions
with downloaded inputs, are preserved. Callers must choose distinct paths
for concurrent calls or separate R processes.

## Consequences

Users can retain a compact record of the actual inputs to a successful call
without extra requests, dependencies or changes to analysis objects.
Hashing adds disk reads only when explicitly enabled. A new output path must
be chosen for each call. Historic headers unavailable in caches stay unknown.

The collector uses package-local state for synchronous nested calls, with
deterministic cleanup. Separate R processes do not share that state.
Transformation provenance, remote freshness checks, attaching metadata to
objects (#224), and receipts for stand-alone lower-level download APIs are
outside this change.

## Alternatives considered

- A new getGEO flag: simpler to discover, but conflicts with ADR-0003.
- Attach metadata or attributes to each return value: changes object contracts
  and overlaps #224.
- Fetch missing dates with additional requests: adds network work, and today's
  remote headers would not establish the origin of an older local file.
- Guess a URL for local files: confuses expected locations with observed origin.
- Write a receipt on every exit: could present a failed parse or partial fetch
  as a successful retrieval.
- Overwrite a previous receipt: risks losing the earlier retrieval record.
