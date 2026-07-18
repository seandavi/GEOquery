# Resolve a GEO accession to its SRA accessions

Maps a GEO Series (`GSE`) or Sample (`GSM`) to the underlying Sequence
Read Archive (SRA) accessions — study (`SRP`), experiment (`SRX`), run
(`SRR`), and sample (`SRS`) — by cross-referencing NCBI Entrez (`gds`
-\> `sra`). This is the accession-linking bridge to raw sequencing data;
it downloads no sequence files.

## Usage

``` r
geoToSRA(geo)
```

## Arguments

- geo:

  character(1), a GEO accession (`GSE...` or `GSM...`).

## Value

A `data.frame` with one row per SRA run and columns `geo_accession`,
`srp`, `srx`, `srr`, `srs`. Returns a zero-row data.frame (same columns)
when the accession has no linked SRA records.

## Details

Uses the NCBI Entrez API through rentrez. An `ENTREZ_KEY` raises the
rate limit; see
[searchGEO](http://seandavi.github.io/GEOquery/reference/searchGEO.md)
for how to set one.

## See also

[searchGEO](http://seandavi.github.io/GEOquery/reference/searchGEO.md)

## Examples

``` r
if (FALSE) { # \dontrun{
geoToSRA("GSE164073")
} # }
```
