# ADR-0008: Bridge GEO accessions to raw sequencing data (SRA linking + ENA files)

- Status: accepted
- Date: 2026-07-18
- Deciders: Sean Davis

## Context

GEOquery retrieves *processed* GEO data and NCBI-*computed* RNA-seq counts, but
offered no path from a GEO accession to the underlying **raw sequencing data**.
Users who needed FASTQ/BAM had to leave R for Python tooling (`pysradb`, `ffq`)
or the abandoned `SRAdb`. Two capabilities are missing and are naturally
separate (issues #219 and #220):

1. **Accession linking** — mapping a GEO Series/Sample (`GSE`/`GSM`) to its SRA
   accessions (study `SRP`, experiment `SRX`, run `SRR`, sample `SRS`). NCBI is
   the authority here: GEO→SRA links live in Entrez.
2. **File resolution** — turning those SRA accessions into concrete download
   URLs with sizes and checksums. The European Nucleotide Archive (ENA) mirrors
   SRA and exposes exactly this via its `filereport` API (FASTQ URLs, md5,
   bytes) far more directly than NCBI does.

These add a *new external API surface* (Entrez `elink`/`efetch`, ENA
`filereport`) and a *new return type* (a tidy accession/file table), which is
why they warrant an ADR even though each is small.

## Decision

We will add a two-layer raw-data bridge, keeping GEOquery a *pointer* to raw
data — it links and resolves URLs but does not host, download, or process
sequence files.

- **Layer 1 — SRA linking (`geoToSRA()`, #219).** Resolve `GSE`/`GSM` → SRA via
  Entrez: `esearch` the `gds` UID, `elink` `gds`→`sra`, then `efetch`
  `rettype=runinfo`. The runinfo CSV is the single richest response — one row
  per run with Study/Experiment/Sample/Run accessions as columns — so we parse
  it into a tidy `data.frame(geo_accession, srp, srx, srr, srs)`. The pure
  parser (`.parse_sra_runinfo()`) is unit-tested offline against a saved
  runinfo fixture; the network wrapper is integration-only.
- **Layer 2 — ENA file resolver (#220).** Consume Layer 1's SRA accessions and
  query ENA `filereport` for FASTQ URLs, md5, and byte sizes. Verified downloads
  reuse the checksum helper from #222 (`.verify_md5()`).

## Consequences

- GEOquery now completes the bridge from GEO to raw data without hosting or
  reprocessing it — a frequently requested capability.
- New dependency on the availability and response shape of two external APIs
  (Entrez, ENA). We isolate the fragile parsing in pure functions with fixture
  tests so format drift surfaces as a failing unit test, and keep live calls in
  the integration suite.
- `geoToSRA()` returns run-level rows (finest granularity); callers aggregate to
  SRP/SRX as needed. Aggregation helpers can follow if demand appears (YAGNI).
- OmicIDX bulk indexing is intentionally out of scope; this uses NCBI/ENA
  directly. An OmicIDX accelerator remains a possible later, optional add.

## Alternatives considered

- **Parse the SRA link out of the GSM SOFT record** (`!Sample_relation = SRA:
  ...`). Works for a lone GSM but not a Series without fetching every sample,
  and gives only the SRX, not the run/sample/study set. Entrez `elink` is the
  authoritative, complete mapping.
- **NCBI-only file resolution** (efetch full SRA XML). ENA `filereport` returns
  URLs + md5 + bytes in one flat TSV; NCBI requires assembling download paths
  and offers no md5. ENA is the lazier, richer source (mirrors SRA).
- **Vendor `pysradb`/`ffq` semantics wholesale.** Out of scope — we only need
  the linking + URL resolution, staying in-lane as a retrieval library.
