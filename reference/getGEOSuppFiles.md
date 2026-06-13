# Get Supplemental Files from GEO

NCBI GEO allows supplemental files to be attached to GEO Series (GSE),
GEO platforms (GPL), and GEO samples (GSM). This function 'knows' how to
get these files based on the GEO accession. No parsing of the downloaded
files is attempted, since the file format is not generally knowable by
the computer.

## Usage

``` r
getGEOSuppFiles(
  GEO,
  makeDirectory = TRUE,
  baseDir = getwd(),
  fetch_files = TRUE,
  filter_regex = NULL,
  quiet = getOption("GEOquery.quiet", FALSE)
)
```

## Arguments

- GEO:

  A GEO accession number such as GPL1073 or GSM1137

- makeDirectory:

  Should a 'subdirectory' for the downloaded files be created? Default
  is TRUE. If FALSE, the files will be downloaded directly into the
  baseDir.

- baseDir:

  The base directory for the downloads. Default is the current working
  directory.

- fetch_files:

  logical(1). If TRUE, then actually download the files. If FALSE, just
  return the filenames that would have been downloaded. Useful for
  testing and getting a list of files without actual download.

- filter_regex:

  A character(1) regular expression that will be used to filter the
  filenames from GEO to limit those files that will be downloaded. This
  is useful to limit to, for example, bed files only.

- quiet:

  logical(1). If TRUE, suppress informational messages such as "No
  supplemental files found" and "Using locally cached version". Defaults
  to the \`GEOquery.quiet\` option, or FALSE.

## Value

If fetch_files=TRUE, a data frame is returned invisibly with rownames
representing the full path of the resulting downloaded files and the
records in the data.frame the output of file.info for each downloaded
file. If fetch_files=FALSE, a data.frame of URLs and filenames is
returned.

## Details

Again, just a note that the files are simply downloaded.

## Examples
