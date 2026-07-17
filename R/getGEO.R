#' Get a GEO object from NCBI or file
#' 
#' This function is the main user-level function in the GEOquery package.  It
#' directs the download (if no filename is specified) and parsing of a GEO SOFT
#' format file into an R data structure specifically designed to make access to
#' each of the important parts of the GEO SOFT format easily accessible.
#' 
#' getGEO functions to download and parse information available from NCBI GEO
#' (\url{http://www.ncbi.nlm.nih.gov/geo}).  Here are some details about what
#' is avaible from GEO.  All entity types are handled by getGEO and essentially
#' any information in the GEO SOFT format is reflected in the resulting data
#' structure.
#' 
#' From the GEO website:
#' 
#' The Gene Expression Omnibus (GEO) from NCBI serves as a public repository
#' for a wide range of high-throughput experimental data. These data include
#' single and dual channel microarray-based experiments measuring mRNA, genomic
#' DNA, and protein abundance, as well as non-array techniques such as serial
#' analysis of gene expression (SAGE), and mass spectrometry proteomic data. At
#' the most basic level of organization of GEO, there are three entity types
#' that may be supplied by users: Platforms, Samples, and Series.
#' Additionally, there is a curated entity called a GEO dataset.
#' 
#' A Platform record describes the list of elements on the array (e.g., cDNAs,
#' oligonucleotide probesets, ORFs, antibodies) or the list of elements that
#' may be detected and quantified in that experiment (e.g., SAGE tags,
#' peptides). Each Platform record is assigned a unique and stable GEO
#' accession number (GPLxxx). A Platform may reference many Samples that have
#' been submitted by multiple submitters.
#' 
#' A Sample record describes the conditions under which an individual Sample
#' was handled, the manipulations it underwent, and the abundance measurement
#' of each element derived from it. Each Sample record is assigned a unique and
#' stable GEO accession number (GSMxxx). A Sample entity must reference only
#' one Platform and may be included in multiple Series.
#' 
#' A Series record defines a set of related Samples considered to be part of a
#' group, how the Samples are related, and if and how they are ordered. A
#' Series provides a focal point and description of the experiment as a whole.
#' Series records may also contain tables describing extracted data, summary
#' conclusions, or analyses. Each Series record is assigned a unique and stable
#' GEO accession number (GSExxx).
#' 
#' GEO DataSets (GDSxxx) are curated sets of GEO Sample data. A GDS record
#' represents a collection of biologically and statistically comparable GEO
#' Samples and forms the basis of GEO's suite of data display and analysis
#' tools. Samples within a GDS refer to the same Platform, that is, they share
#' a common set of probe elements. Value measurements for each Sample within a
#' GDS are assumed to be calculated in an equivalent manner, that is,
#' considerations such as background processing and normalization are
#' consistent across the dataset. Information reflecting experimental design is
#' provided through GDS subsets.
#' 
#' @param GEO A character string representing a GEO object for download and
#' parsing.  (eg., 'GDS505','GSE2','GSM2','GPL96')
#' @param filename The filename of a previously downloaded GEO SOFT format file
#' or its gzipped representation (in which case the filename must end in .gz).
#' Either one of GEO or filename may be specified, not both.  GEO series matrix
#' files are also handled.  Note that since a single file is being parsed, the
#' return value is not a list of esets, but a single eset when GSE matrix files
#' are parsed.
#' @param destdir The destination directory for any downloads.  Defaults to the
#' architecture-dependent tempdir.  You may want to specify a different
#' directory if you want to save the file for later use.  Doing so is a good
#' idea if you have a slow connection, as some of the GEO files are HUGE!
#' @param GSElimits This argument can be used to load only a contiguous subset
#' of the GSMs from a GSE.  It should be specified as a vector of length 2
#' specifying the start and end (inclusive) GSMs to load.  This could be useful
#' for splitting up large GSEs into more manageable parts, for example.
#' @param GSEMatrix A boolean telling GEOquery whether or not to use GSE Series
#' Matrix files from GEO.  The parsing of these files can be many
#' orders-of-magnitude faster than parsing the GSE SOFT format files.  Defaults
#' to TRUE, meaning that the SOFT format parsing will not occur; set to FALSE
#' if you for some reason need other columns from the GSE records.
#' @param AnnotGPL A boolean defaulting to FALSE as to whether or not to use
#' the Annotation GPL information.  These files are nice to use because they
#' contain up-to-date information remapped from Entrez Gene on a regular basis.
#' However, they do not exist for all GPLs; in general, they are only available
#' for GPLs referenced by a GDS
#' @param getGPL A boolean defaulting to TRUE as to whether or not to download
#' and include GPL information when getting a GSEMatrix file.  You may want to
#' set this to FALSE if you know that you are going to annotate your
#' featureData using Bioconductor tools rather than relying on information
#' provided through NCBI GEO.  Download times can also be greatly reduced by
#' specifying FALSE.
#' @param parseCharacteristics A boolean defaulting to TRUE as to whether or not
#' to parse the characteristics information (if available) for a GSE Matrix file.
#' Set this to FALSE if you experience trouble while parsing the characteristics.
#' @param returnType One of "SummarizedExperiment" (default) or
#' "ExpressionSet". For GSE Series Matrix results, controls whether each entity
#' is returned as a SummarizedExperiment or an ExpressionSet. SOFT-format
#' results (GDS/GPL/GSM/GSE S4 objects) are unaffected. As of this release the
#' default is "SummarizedExperiment"; pass returnType = "ExpressionSet" for the
#' previous behavior.
#' @param encoding Optional character encoding for reading the downloaded GEO
#' file, one of "unknown" (the default; let the reader auto-detect), "UTF-8", or
#' "Latin-1". Most GEO files are UTF-8/ASCII, but a few carry Latin-1 bytes that
#' are otherwise mis-decoded; set \code{encoding = "Latin-1"} for those. Applies
#' for the duration of the call only. Equivalent to setting
#' \code{options(GEOquery.encoding = ...)} globally.
#' @param token Optional NCBI GEO reviewer access token (character(1)) for
#' fetching a private/embargoed record. Obtain it from the "Reviewer access"
#' link on the private GSE's GEO page. Because private records are not published
#' to the GEO FTP tree, supplying a token forces the SOFT (\code{acc.cgi}) path:
#' for a GSE this returns a \code{GSE} S4 object (as with
#' \code{GSEMatrix = FALSE}), not a \code{SummarizedExperiment} /
#' \code{ExpressionSet}. See ADR-0007.
#' @return An object of the appropriate class (GDS, GPL, GSM, or GSE) is
#' returned.  If the GSEMatrix option is used, then a list of
#' SummarizedExperiment objects is returned by default (or ExpressionSet
#' objects if \code{returnType = "ExpressionSet"}), one for each SeriesMatrix
#' file associated with the GSE accession.
#' @section Warning : Some of the files that are downloaded, particularly those
#' associated with GSE entries from GEO are absolutely ENORMOUS and parsing
#' them can take quite some time and memory.  So, particularly when working
#' with large GSE entries, expect that you may need a good chunk of memory and
#' that coffee may be involved when parsing....
#' 
#' @importFrom readr problems
#' 
#' @author Sean Davis
#' @seealso \code{\link{getGEOfile}}
#' @keywords IO
#' @examples
#' \dontrun{
#' 
#' gds <- getGEO('GDS10')
#' gds
#'
#' gse <- getGEO('GSE10')
#' # Returns a list, so look at first item
#' 
#' gse[[1]]
#' 
#' }
#' @export
getGEO <- function(GEO = NULL, filename = NULL, destdir = tempdir(), GSElimits = NULL,
    GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = TRUE, parseCharacteristics = TRUE,
    returnType = c("SummarizedExperiment", "ExpressionSet"), encoding = NULL,
    token = NULL) {
    returnType_default <- missing(returnType)
    returnType <- match.arg(returnType)
    if (!is.null(token) && (!is.character(token) || length(token) != 1L)) {
        stop("'token' must be a single character string (an NCBI GEO reviewer access token)")
    }
    # Optional per-call character encoding for the underlying GEO file reads
    # (data.table::fread). Sets the GEOquery.encoding option for the duration of
    # this call only, restoring the previous value on exit. Useful for the
    # occasional non-UTF-8 GEO record (#148).
    if (!is.null(encoding)) {
        encoding <- match.arg(encoding, c("unknown", "UTF-8", "Latin-1"))
        old_encoding <- getOption("GEOquery.encoding")
        options(GEOquery.encoding = encoding)
        on.exit(options(GEOquery.encoding = old_encoding), add = TRUE)
    }
    con <- NULL
    if (!is.null(GSElimits)) {
        if (length(GSElimits) != 2) {
            stop("GSElimits should be an integer vector of length 2, like (1,10) to include GSMs 1 through 10")
        }
    }
    if (is.null(GEO) & is.null(filename)) {
        stop("You must supply either a filename of a GEO file or a GEO accession")
    }
    if (is.null(filename)) {
        GEO <- toupper(GEO)
        geotype <- toupper(substr(GEO, 1, 3))
        # A reviewer token means a private record, which has no Series Matrix on
        # the FTP tree; fall through to the SOFT (acc.cgi) path, which returns a
        # GSE S4 object rather than an ExpressionSet/SummarizedExperiment (#154).
        if (GSEMatrix & geotype == "GSE" & is.null(token)) {
            ret <- getAndParseGSEMatrices(GEO, destdir, AnnotGPL = AnnotGPL, getGPL = getGPL,
                parseCharacteristics = parseCharacteristics)
            return(.applyReturnType(ret, returnType, returnType_default))
        }
        filename <- getGEOfile(GEO, destdir = destdir, AnnotGPL = AnnotGPL, token = token)
    }
    ret <- parseGEO(filename, GSElimits, destdir, AnnotGPL = AnnotGPL, getGPL = getGPL,
        parseCharacteristics = parseCharacteristics)
    return(.applyReturnType(ret, returnType, returnType_default))
}

# Apply the requested return type to a parsed result. ExpressionSet results
# (the GSE Series Matrix path) may be coerced to SummarizedExperiment; SOFT S4
# objects (GSE/GSM/GPL/GDS) are returned unchanged. See ADR-0002 (#168).
.applyReturnType <- function(ret, returnType, notify_default = FALSE) {
    is_eset <- function(x) methods::is(x, "ExpressionSet")
    contains_eset <- (is.list(ret) && length(ret) > 0 && is_eset(ret[[1]])) || is_eset(ret)

    if (notify_default && returnType == "SummarizedExperiment" && contains_eset) {
        rlang::inform(
            paste0(
                "getGEO() now returns SummarizedExperiment objects by default. ",
                "Pass returnType = 'ExpressionSet' for the previous behavior."
            ),
            .frequency = "once", .frequency_id = "geoquery_returnType_default"
        )
    }

    if (returnType == "ExpressionSet") {
        return(ret)
    }
    coerce_one <- function(x) if (is_eset(x)) as_SummarizedExperiment(x) else x
    if (is.list(ret)) {
        return(lapply(ret, coerce_one))
    }
    coerce_one(ret)
}

#' Coerce a GEOquery ExpressionSet to a SummarizedExperiment
#'
#' A thin wrapper around
#' \code{SummarizedExperiment::makeSummarizedExperimentFromExpressionSet()} used
#' by \code{getGEO(..., returnType = "SummarizedExperiment")}, and available
#' directly so existing ExpressionSet results can be modernized without
#' re-downloading.
#'
#' @param eset An \code{ExpressionSet}, e.g. an element returned by
#'   \code{getGEO()} for a GSE Series Matrix file.
#' @return A \code{SummarizedExperiment}.
#' @examples
#' \dontrun{
#'   gse <- getGEO("GSE2553")[[1]]
#'   se <- as_SummarizedExperiment(gse)
#' }
#' @export
as_SummarizedExperiment <- function(eset) {
    if (!methods::is(eset, "ExpressionSet")) {
        stop("'eset' must be an ExpressionSet")
    }
    SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(eset)
}
