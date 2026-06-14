# Single-cell example datasets (GEO)

Real GEO accessions used to design and test the single-cell readers
(`geoSingleCellManifest()`, `geoSingleCellUnits()`, `getGEOSingleCell()`,
`readGEOSingleCell()`; see ADR-0004). The goal is to exercise **all three
import pathways** and the **layout variants** (per-sample vs whole-study,
series-level vs `_RAW.tar` vs per-GSM suppl dir, single vs multi-platform)
without over- or under-designing the API.

Maintainer reference only — this directory is excluded from the package build
(`.Rbuildignore`). All facts below were confirmed against live NCBI GEO
(verified 2026-06-14); GEO layouts can change.

## Import pathways

| Pathway | Format label | Reader | Backend |
|---|---|---|---|
| 10x Matrix Market triplet | `10x_mtx` | `readGEOSingleCell()` | TENxIO |
| 10x HDF5 (CellRanger `.h5`) | `10x_h5` | `readGEOSingleCell()` | TENxIO |
| AnnData `.h5ad` | `h5ad` | `readGEOSingleCell()` | anndataR |
| loom | `loom` | — | out of scope (native pkg) |
| Seurat `.rds` | `rds` | — | out of scope (native pkg) |

## Example accessions

| Accession | Format | Pathway | Layout | Where files live | Platforms | What it exercises |
|---|---|---|---|---|---|---|
| **GSE132771** | `10x_mtx` | TENxIO | per-sample | `_RAW.tar` **and** per-GSM suppl | **2** — GPL21103 (mouse, 8), GPL24676 (human, 16) | GSE→GSM fallback; `by="platform"` grouping; cross-platform combine error; mixed CellRanger v2 `genes.tsv` / v3 `features.tsv` → differing rowData |
| **GSE145926** | `10x_h5` | TENxIO | per-sample | `_RAW.tar` **and** per-GSM suppl (`*_filtered_feature_bc_matrix.h5`) | 1 (12) | clean per-sample 10x HDF5. **Reader verified** (GSM4339769 → SCE 33539×6249) |
| **GSE122960** | `10x_h5` | TENxIO | per-sample | `_RAW.tar` (`*_filtered_gene_bc_matrices_h5.h5`) | 1 (17) | second 10x HDF5 example; alt filename |
| **GSE310450** | `h5ad` | anndataR | **whole-study** | series-level (`GSE310450_WTS_LACZ.h5ad`, no GSM) | 1 | ~62 MB modern anndata, genuinely single-cell. **Reader verified** (→ SCE 14021×8603) |
| **GSE328481** | `h5ad` | anndataR | **per-sample** | per-GSM suppl (`GSM9684206_*.h5ad`, one per sample) | 1 | per-sample h5ad layout (files ~1 GB each, so large to read) |
| GSE312831 | `h5ad` | anndataR | whole-study | series-level (`GSE312831_*.h5ad`) | 1 | tiny (~2 MB) but **bulk** (`bulkseq`, 6 columns) — readable, not single-cell; not used as the reader example |
| **GSE161228** | `h5ad` | anndataR | **whole-study** | series-level loose (`GSE161228_*.h5ad.gz`, no GSM) | 1 | already-combined single files; `sample = NA`. **Classification + gunzip verified, reader fails**: files are anndata<0.8.0, which anndataR cannot import (a 2020 study) |
| **GSE154567** | `10x_h5` + `rds` | TENxIO (h5) | per-sample | `_RAW.tar` (9 `.h5` + 9 `.rds`) | 1 (9) | mixed formats in one study; format selection/preference |
| **GSE150728** | `rds` (Seurat) | — | whole-study | series-level + `_RAW.tar` | 1 (13) | manifest classification of out-of-scope format |
| **GSE131907** | `rds` (Seurat) | — | whole-study | series-level (2 `.rds`) | 1 (58) | manifest classification of out-of-scope format |

## Layout dimensions (why these were chosen)

- **Per-sample vs whole-study.** Per-sample (GSE132771, GSE145926) → one object
  per GSM, grouping/combine is meaningful. Whole-study (GSE161228 h5ad,
  GSE131907 rds) → one file already holds everything, so `by`/combine is a
  no-op and `sample` is `NA`.
- **File location.** Loose at the series level (GSE161228), inside `_RAW.tar`
  (most), or in per-GSM suppl dirs (GSE132771, GSE145926). The manifest's
  GSE→GSM fallback exists for the last case.
- **Single vs multi-platform.** GSE132771 is the key multi-platform case (mouse
  + human); it's why `by = "platform"` is the right combine boundary rather than
  one object for the whole series.
- **Annotation heterogeneity.** Within GSE132771's human platform, samples mix
  CellRanger v2 (`genes.tsv` → rowData `ID, Symbol`) and v3 (`features.tsv` →
  `ID, Symbol, Type`), so `combine`/`by="all"` must reconcile rowData columns.

## Reader-pathway coverage (verified live 2026-06-14)

All three import pathways verified end-to-end:

- `10x_mtx` (TENxIO) — **verified** (GSE132771 GSM3891612 → SCE 27998×4245).
- `10x_h5` (TENxIO) — **verified** (GSE145926 GSM4339769 → SCE 33539×6249).
- `h5ad` (anndataR) — **verified** (GSE310450 → SCE 14021×8603, real
  single-cell). Older h5ad (anndata&lt;0.8.0, e.g. GSE161228) cannot be read;
  prefer recent submissions.

## Finding more examples (OmicIDX GEO parquet)

GEO metadata, including each accession's supplemental file list, is published as
DuckDB-queryable parquet — far easier than scraping FTP:

- `https://data-omicidx.cancerdatasci.org/geo/parquet/geo_series.parquet`
- `https://data-omicidx.cancerdatasci.org/geo/parquet/geo_samples.parquet`

`supplemental_files` is a `varchar[]` of FTP URLs. To find recent readable h5ad:

```sql
INSTALL httpfs; LOAD httpfs;
SELECT accession, submission_date,
       list_filter(supplemental_files, x -> x ILIKE '%h5ad%') AS h5ad
FROM read_parquet('https://data-omicidx.cancerdatasci.org/geo/parquet/geo_series.parquet')
WHERE array_to_string(supplemental_files, ';') ILIKE '%.h5ad%'
ORDER BY submission_date DESC LIMIT 25;
```

`geo_series` carries `sample_id[]` (its GSMs), `sample_organism[]`, `platform_id[]`;
`geo_samples` is keyed by GSM. The parquet has no file sizes — `curl -I` the URL.

## Gaps still wanted

- A **loom** study — manifest classification only (reader out of scope).
- A confirmed **multi-platform h5/h5ad** study — to test `by="platform"` beyond
  the 10x_mtx case (all h5/h5ad examples so far are single-platform).
