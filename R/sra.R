# Resolve GEO accessions to the underlying SRA accessions (#219, ADR-0008).
#
# GEOquery retrieves processed data and NCBI-computed RNA-seq counts, but had no
# way to get from a GEO accession to the raw-sequencing SRA accessions (study
# SRP, experiment SRX, run SRR, sample SRS). This is the *linking* layer only:
# it cross-references accessions via NCBI Entrez and returns a tidy table. It
# does not download sequencing data (that is the ENA raw-file resolver, #220,
# which consumes this table).
#
# Pipeline: accession -> gds UID (esearch) -> sra UIDs (elink gds->sra) ->
# runinfo CSV (efetch db=sra rettype=runinfo). The runinfo CSV is the single
# richest Entrez response: one row per run, with the Study/Experiment/Sample/Run
# accessions all present as columns. `.parse_sra_runinfo()` (the offline-tested
# core) maps those columns to the tidy `{geo_accession, srp, srx, srr, srs}`.

#' Resolve a GEO accession to its SRA accessions
#'
#' Maps a GEO Series (`GSE`) or Sample (`GSM`) to the underlying Sequence Read
#' Archive (SRA) accessions — study (`SRP`), experiment (`SRX`), run (`SRR`),
#' and sample (`SRS`) — by cross-referencing NCBI Entrez (`gds` -> `sra`). This
#' is the accession-linking bridge to raw sequencing data; it downloads no
#' sequence files.
#'
#' Uses the NCBI Entrez API through \pkg{rentrez}. An `ENTREZ_KEY` raises the
#' rate limit; see [searchGEO] for how to set one.
#'
#' @param geo character(1), a GEO accession (`GSE...` or `GSM...`).
#'
#' @return A `data.frame` with one row per SRA run and columns
#'   `geo_accession`, `srp`, `srx`, `srr`, `srs`. Returns a zero-row data.frame
#'   (same columns) when the accession has no linked SRA records.
#'
#' @seealso [searchGEO]
#'
#' @examples
#' \dontrun{
#' geoToSRA("GSE164073")
#' }
#'
#' @importFrom rentrez entrez_search entrez_link entrez_fetch
#' @export
geoToSRA <- function(geo) {
  if (length(geo) != 1L || !is.character(geo)) {
    .abort_bad_accession(geo)
  }
  geo <- toupper(geo)
  if (!grepl("^GS[EM]\\d+$", geo)) {
    .abort_bad_accession(geo)
  }
  uids <- rentrez::entrez_search("gds", paste0(geo, "[ACCN]"), retmax = 50L)$ids
  if (!length(uids)) {
    return(.empty_sra_table())
  }
  linked <- rentrez::entrez_link(dbfrom = "gds", db = "sra", id = uids)
  sra_ids <- linked$links$gds_sra
  if (!length(sra_ids)) {
    return(.empty_sra_table())
  }
  runinfo <- rentrez::entrez_fetch(
    db = "sra", id = sra_ids,
    rettype = "runinfo", retmode = "text"
  )
  .parse_sra_runinfo(runinfo, geo)
}

# Column names in an SRA "runinfo" CSV that carry each accession level.
.SRA_RUNINFO_COLS <- c(srp = "SRAStudy", srx = "Experiment", srr = "Run", srs = "Sample")

# Parse an SRA runinfo CSV (text) into the tidy {geo_accession, srp, srx, srr,
# srs} table. Offline-testable core of geoToSRA(); takes no network.
.parse_sra_runinfo <- function(runinfo, geo_accession) {
  dt <- data.table::fread(text = runinfo, sep = ",", header = TRUE, fill = TRUE)
  # An empty runinfo body (header only, or nothing) -> no runs.
  if (!nrow(dt) || !all(.SRA_RUNINFO_COLS %in% names(dt))) {
    return(.empty_sra_table())
  }
  out <- data.frame(
    geo_accession = geo_accession,
    srp = as.character(dt[[.SRA_RUNINFO_COLS[["srp"]]]]),
    srx = as.character(dt[[.SRA_RUNINFO_COLS[["srx"]]]]),
    srr = as.character(dt[[.SRA_RUNINFO_COLS[["srr"]]]]),
    srs = as.character(dt[[.SRA_RUNINFO_COLS[["srs"]]]]),
    stringsAsFactors = FALSE
  )
  # Drop blank runs (trailing empty CSV lines parse to "").
  out[nzchar(out$srr), , drop = FALSE]
}

.empty_sra_table <- function() {
  data.frame(
    geo_accession = character(0), srp = character(0), srx = character(0),
    srr = character(0), srs = character(0), stringsAsFactors = FALSE
  )
}
