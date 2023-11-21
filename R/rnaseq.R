.gse_download_types <- readr::read_csv(
  I(
    "id,description
    soft,SOFT formatted metadata and data
    matrix,GSEmatrix formatted data
    miniml,MINiML formatted metadata and data
    suppl_archive,Supplementary files as tar archive
    raw_counts,Raw RNA-seq counts
    norm_counts_fpkm,Normalized RNA-seq counts (FPKM)
    norm_counts_tpm,Normalized RNA-seq counts (TPM)"
  ),
  show_col_types = FALSE,
  progress = FALSE
)

#' Get available files for a GEO accession number
#'
#' @param gse GEO accession number
#'
#' @return A list of data frames
#'
#' @examples
#'
#' available_gse_files("GSE83322")
#'
#' @export
available_gse_files <- function(gse) {
  # TODO: check if gse is valid
  # TODO: need to add support for downloading annotation for RNA-seq
  # TODO: should add file size and date where available

  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gse)
  page <- rvest::read_html(url)
  anchors <- page |>
    rvest::html_nodes(xpath = "//li/a")
  links <- anchors |>
    rvest::html_attr("href")
  links <- sub("^ftp", "https", links)
  link_ids <- anchors |>
    rvest::html_attr("id")

  links[!grepl("^https://ftp.ncbi.nlm.nih.gov", links)] <- paste0(
    "https://www.ncbi.nlm.nih.gov",
    links[!grepl("^https://ftp.ncbi.nlm.nih.gov/", links)]
  )

  data.frame(link = links, id = link_ids) |>
    dplyr::filter(grepl("^download", id)) |> # nolint: object_usage_linter.
    dplyr::mutate(id = sub("download_", "", id)) |>
    dplyr::left_join(.gse_download_types, by = "id")
}
