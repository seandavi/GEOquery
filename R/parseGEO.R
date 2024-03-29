.na_strings = c("NA", "null", "NULL", "Null")


#' @importFrom data.table fread
.read_lines <- function(...) {
    data.table::fread(..., sep = "", header = FALSE)[[1]]
}

#' Parse GEO text
#' 
#' Workhorse GEO parsers.
#' 
#' These are probably not useful to the end-user.  Use getGEO to access these
#' functions.  parseGEO simply delegates to the appropriate specific parser.
#' There should be no reason to use the parseGPL, parseGDS, parseGSE, or
#' parseGSM functions directly.
#' 
#' @aliases parseGEO parseGPL parseGSE parseGDS parseGSM
#' @param fname The filename of a SOFT format file.  If the filename ends in
#' .gz, a gzfile() connection is used to read the file directly.
#' @param GSElimits Used to limit the number of GSMs parsed into the GSE
#' object; useful for memory management for large GSEs.
#' @param destdir The destination directory into which files will be saved (to
#' be used for caching)
#' @param AnnotGPL Fetch the annotation GPL if available
#' @param getGPL Fetch the GPL associated with a GSEMatrix entity (should
#' remain TRUE for all normal use cases)
#' @return parseGEO returns an object of the associated type.  For example, if
#' it is passed the text from a GDS entry, a GDS object is returned.
#' @author Sean Davis
#' @seealso \code{\link{getGEO}}
#' @keywords IO
#' 
#' @export
parseGEO <- function(fname, GSElimits, destdir = tempdir(), AnnotGPL = FALSE, getGPL = TRUE) {
    con <- fileOpen(fname)
    first.entity <- findFirstEntity(con)
    close(con)
    ret <- switch(as.character(first.entity[1]), sample = {
        parseGSM(fname)
    }, series = parseGSE(fname, GSElimits), dataset = {
        parseGDS(fname)
    }, platform = {
        parseGPL(fname)
    }, `0` = {
        parseGSEMatrix(fname, destdir = destdir, AnnotGPL = AnnotGPL, getGPL = getGPL)$eset
    }, )
    return(ret)
}

parseGeoMeta <- function(txt) {
    ## leader <- strsplit(grep('!\\w*_',txt,perl=TRUE,value=TRUE)[1],'_')[[1]][1]
    ## pull out only lines that are in the header
    tmp <- txt[grep("!\\w*?_", txt)]
    tmp <- gsub("!\\w*?_", "", tmp)
    first.eq <- regexpr(" = ", tmp)
    tmp <- cbind(substring(tmp, first = 1, last = first.eq - 1), substring(tmp, first = first.eq +
        3))
    tmp <- tmp[tmp[, 1] != "", ]
    header <- split(tmp[, 2], tmp[, 1])
    return(header)
}

parseGDSSubsets <- function(txt) {
    ## takes GDS text as input returns a data frame suitable for inclusion as a
    ## 'Column' slot in a GDS GeoDataTable object
    numSubsets <- length(grep("^\\^subset", txt, ignore.case = TRUE))
    subset.lists <- list()
    if (numSubsets > 0) {
        subset.types <- do.call("rbind", strsplit(txt[grep("!subset_type", txt)], " = "))[,
            2]
        subset.descriptions <- do.call("rbind", strsplit(txt[grep("!subset_description",
            txt)], " = "))[, 2]
        subset.samples <- strsplit(do.call("rbind", strsplit(txt[grep("!subset_sample_id",
            txt)], " = "))[, 2], ",")
        for (i in 1:length(subset.types)) {
            for (j in 1:length(subset.samples[[i]])) {
                subset.lists[[subset.types[i]]][[subset.samples[[i]][j]]] <- subset.descriptions[i]
            }
        }
    }
    sample.descriptions <- do.call("rbind", strsplit(txt[grep("#GSM", txt)], " = "))
    sample.descriptions[, 1] <- gsub("#", "", sample.descriptions[, 1])
    samp.table <- data.frame(sample = as.character(sample.descriptions[, 1]))
    if (length(subset.lists) > 0) {
        for (i in 1:length(subset.lists)) {
            samp.table[match(names(subset.lists[[i]]), samp.table[, 1]), names(subset.lists)[i]] <- as.factor(unlist(subset.lists[[i]]))
        }
        colnames(samp.table) <- c("sample", gsub(" ", ".", names(subset.lists)))
    } else {
        colnames(samp.table) <- "sample"
    }
    samp.table[match(sample.descriptions[, 1], samp.table[, 1]), "description"] <- sample.descriptions[,
        2]
    return(as.data.frame(samp.table))
}


splitOnFirst <- function(x, pattern) {
    patlen <- nchar(pattern)
    matches <- regexpr(pattern, x)
    leftside <- substr(x, start = 1, stop = matches - 1)
    rightside <- substr(x, start = matches + patlen, stop = 1e+07)
    return(data.frame(leftside, rightside))
}


parseGeoColumns <- function(txt) {
    cols <- as.data.frame(splitOnFirst(txt[grep("^#", txt, perl = TRUE)], " = "))
    cols[, 1] <- sub("#", "", as.character(cols[, 1]))
    colnames(cols) <- c("Column", "Description")
    return(cols)
}

### Parse a GSM; the limit is established by the parameter n, used by getGSE to
### limit the number of lines read to only the size of ONE GSM, since a single
### GSE contains many GSMs
#' @importFrom readr read_lines
.parseGSMWithLimits <- function(con, n = NULL) {
    txt <- vector("character")
    i <- 0
    hasDataTable = FALSE
    while (i <- i + 1) {
        tmp <- try(read_lines(con, 1))
        if (inherits(tmp, "try-error") | length(tmp) == 0) {
            i <- i - 1  # counted an extra line
            hasDataTable = FALSE
            break
        }
        txt[i] <- tmp
        if (i == length(txt))
            txt <- c(txt, character(i))  # double vector size
        if (length(grep("!\\w+_table_begin", txt[i], perl = TRUE)) > 0) {
            hasDataTable = TRUE
            break
        }
    }
    txt <- txt[1:i]
    cols <- parseGeoColumns(txt)
    meta <- parseGeoMeta(txt)
    geoDataTable <- new("GEODataTable", columns = data.frame(), table = data.frame())
    if (hasDataTable) {
        nLinesToRead <- NULL
        if (!is.null(n)) {
            nLinesToRead <- n - length(txt)
        }
        dat3 <- fastTabRead(con, n = nLinesToRead)
        geoDataTable <- new("GEODataTable", columns = cols, table = dat3[1:(nrow(dat3) -
            1), ])
    }
    gsm <- new("GSM", header = meta, dataTable = geoDataTable)
}


### This function does a grep on a file by doing a readline in chunks of size
### chunksize.  Return value is a data.frame with the line number of each match
### and the line itself.
filegrep <- function(con, regex, chunksize = 10000) {
    i <- 0
    ret <- NULL
    while (TRUE) {
        lines <- readLines(con, n = chunksize)
        if (length(lines) == 0) {
            break
        }
        foundLines <- grep(regex, lines)
        foundTypes <- lines[foundLines]
        if (length(foundLines) > 0) {
            foundLines <- foundLines + i
            tmp <- data.frame(foundLines = foundLines, foundTypes = foundTypes)
            if (is.null(ret)) {
                ret <- tmp
            } else {
                ret <- rbind(ret, tmp)
            }
        }
        i <- i + length(lines)
    }
    return(ret)
}

#' @importFrom readr read_lines
parseGSE <- function(fname, GSElimits = NULL) {
    gsmlist <- list()
    gpllist <- list()
    GSMcount <- 0
    message("Reading file....")
    lines = .read_lines(fname)
    message("Parsing....")
    lineCounts <- grep("\\^(SAMPLE|PLATFORM)", lines)
    message(sprintf("Found %d entities...", length(lineCounts)))
    ## I close and reopen the file because on Windows, the seek function is
    ## pretty much guaranteed to NOT work This gets the header information for
    ## the GSE
    a <- lines[1:(lineCounts[1] - 1)]
    header = parseGeoMeta(a)
    lineCounts = c(lineCounts, length(lines))
    ## parse the actual entities, now
    for (j in 1:(length(lineCounts) - 1)) {
        tmp <- strsplit(lines[lineCounts[j]], " = ")[[1]]
        accession <- tmp[2]
        message(sprintf("%s (%d of %d entities)", accession, j, length(lineCounts)))
        entityType <- tolower(sub("\\^", "", tmp[1]))
        if (entityType == "sample") {
            gsmlist[[accession]] <- .parseGSMTxt(lines[lineCounts[j]:(lineCounts[j +
                1] - 1)])
        }
        if (entityType == "platform") {
            gpllist[[accession]] <- .parseGPLTxt(lines[lineCounts[j]:(lineCounts[j +
                1] - 1)])
        }
    }
    return(new("GSE", header = header, gsms = gsmlist, gpls = gpllist))
}


findFirstEntity <- function(con) {
    while (TRUE) {
        line <- suppressWarnings(readLines(con, 100))
        if (length(line) == 0)
            return(0)
        entity.line <- grep("^\\^(DATASET|SAMPLE|SERIES|PLATFORM|ANNOTATION)", line,
            ignore.case = TRUE, value = TRUE, perl = TRUE)
        entity.line <- gsub("annotation", "platform", entity.line, ignore.case = TRUE)
        if (length(entity.line) > 0) {
            ret <- c(tolower(sub("\\^", "", strsplit(entity.line, " = ")[[1]][1])),
                strsplit(entity.line, " = ")[[1]][2])
            return(ret)
        }
        # catch GSE SeriesMatrix files
        checkline = grep("^!Series_title", line, ignore.case = TRUE, value = TRUE,
            perl = TRUE)
        # Messy, but GSEMatrix files use tab-separation rather than '=' for
        # separation, so if the ' = ' is not present, we presume to have a
        # GSEMatrix file (return 0)
        if (length(checkline) > 0) {
            if (!grepl(" = ", checkline))
                return(0)
        }
    }
    return(0)
}

fastTabRead <- function(con, sep = "\t", header = TRUE, sampleRows = 100, colClasses = NULL,
    n = NULL, quote = "\"", ...) {
    ### Need to read tables quickly, so guess the colclasses on the fly.  This is
    ### a bit dangerous since the first rows might not be representative of the
    ### entire set, but it is SO MUCH FASTER than the alternative, I have to do
    ### it.
    dat3 <- data.frame(NULL)
    numberOfLines <- -1
    if (!is.null(n)) {
        numberOfLines <- n - sampleRows
    }
    if (is.null(colClasses)) {
        if (!is.null(n)) {
            sampleRows <- min(sampleRows, n)
        }
        dat1 <- read.table(con, sep = sep, header = header, nrows = sampleRows, fill = TRUE,
            check.names = FALSE, quote = quote, comment.char = "", na.strings = c("NA",
                "null", "NULL", "Null"), ...)
        colclasses <- apply(dat1, 2, class)
        colclasses[1] <- "factor"
        dat2 <- try(read.delim(con, sep = sep, colClasses = colclasses, header = FALSE,
            quote = quote, comment.char = "", na.strings = c("NA", "null", "NULL",
                "Null"), nrows = numberOfLines, ...), silent = TRUE)
        if (inherits(dat2, "try-error")) {
            dat3 = dat1
        } else {
            colnames(dat2) <- colnames(dat1)
            dat3 <- rbind(dat1, dat2)
        }
    } else {
        dat3 <- read.delim(con, sep = sep, colClasses = colClasses, header = header,
            quote = quote, comment.char = "", na.strings = c("NA", "null", "NULL",
                "Null"), nrows = numberOfLines, ...)
    }
    return(dat3)
}

.genericGEOTableParser <- function(txt) {
    # Find first line of the table
    has_table = TRUE
    tbl_begin = grep("!\\w+_table_begin", txt, perl = TRUE)[1] + 1
    # No data table at all, so everything is header info
    if (is.na(tbl_begin)) {
        tbl_begin = length(txt)
        has_table = FALSE
    }

    if (has_table) {
        metadata_rows <- seq_len(tbl_begin - 2)
    } else {
        metadata_rows <- seq_along(txt)
    }


    # Find last line of the table Edge case is that the record does not have a
    # table end line, in which case we use all lines
    tbl_end = grep("!\\w+_table_end", txt, perl = TRUE)
    if (length(tbl_end) == 0) {
        tbl_end = length(txt)
    } else {
        tbl_end = tbl_end[1] - 1
    }
    dat3 = data.frame()
    if (has_table) {
        dat3 <- data.table::fread(text = txt[tbl_begin:tbl_end], na = .na_strings)
    }
    return(list(meta_text = txt[metadata_rows], data_frame = dat3))
}

#' @importFrom data.table fread
parseGDS <- function(fname) {
    txt = data.table::fread(fname, sep = "")[[1]]

    parser_results = .genericGEOTableParser(txt)

    cols <- parseGDSSubsets(parser_results$meta_text)
    meta <- parseGeoMeta(parser_results$meta_text)
    geoDataTable <- new("GEODataTable", columns = cols, table = parser_results$data_frame)
    gds <- new("GDS", header = meta, dataTable = geoDataTable)
    return(gds)
}

.parseGPLWithLimits <- function(con, n = NULL) {
    txt <- vector("character")
    i <- 0
    hasDataTable = FALSE
    while (i <- i + 1) {
        tmp <- try(readLines(con, 1))
        if (inherits(tmp, "try-error") || length(tmp) == 0 || (!is.null(n) && i ==
            n)) {
            i <- i - 1  # counted an extra line
            hasDataTable = FALSE
            break
        }
        txt[i] <- tmp
        if (i == length(txt))
            txt <- c(txt, character(i))  # double vector size
        if (length(grep("!\\w+_table_begin", txt[i], perl = TRUE)) > 0) {
            hasDataTable = TRUE
            break
        }
    }
    txt <- txt[1:i]
    cols <- parseGeoColumns(txt)
    meta <- parseGeoMeta(txt)
    geoDataTable <- new("GEODataTable", columns = data.frame(), table = data.frame())
    if (hasDataTable) {
        nLinesToRead <- NULL
        if (!is.null(n)) {
            nLinesToRead <- n - length(txt)
        }
        dat3 <- fastTabRead(con, n = nLinesToRead, quote = "")
        geoDataTable <- new("GEODataTable", columns = cols, table = dat3[1:(nrow(dat3) -
            1), ])
    }
    gpl <- new("GPL", header = meta, dataTable = geoDataTable)
}


#' @importFrom readr read_tsv
.parseGSMTxt <- function(txt) {
    parser_results = .genericGEOTableParser(txt)

    cols <- parseGeoColumns(parser_results$meta_text)
    meta <- parseGeoMeta(parser_results$meta_text)
    geoDataTable <- new("GEODataTable", columns = cols, table = parser_results$data_frame)
    geo <- new("GSM", header = meta, dataTable = geoDataTable)
    return(geo)
}


#' @importFrom data.table fread
parseGSM <- function(fname) {
    txt = data.table::fread(fname, sep = "")[[1]]
    # read_lines reads separate blank lines on windows, so remove them before
    # proceeding. NOT doing so results in the Table header being repeated as the
    # first line of the Table and test failures galore.
    txt = txt[txt != ""]
    return(.parseGSMTxt(txt))
}

### In memory cache for GPL objects parsed from locally cached versions of GPL
### SOFT files.  It is disabled by default with
### options('GEOquery.inmemory.gpl'=FALSE).
GPLcache <- new.env(parent = emptyenv())


#' @importFrom readr read_tsv
.parseGPLTxt <- function(txt) {

    parser_results = .genericGEOTableParser(txt)

    cols <- parseGeoColumns(parser_results$meta_text)
    meta <- parseGeoMeta(parser_results$meta_text)
    geoDataTable <- new("GEODataTable", columns = cols, table = parser_results$data_frame)
    geo <- new("GPL", header = meta, dataTable = geoDataTable)
    return(geo)
}


#' @importFrom readr read_lines
parseGPL <- function(fname) {
    if (getOption("GEOquery.inmemory.gpl")) {
        info <- file.info(fname, extra_cols = FALSE)
        cache <- get0(fname, envir = GPLcache, inherits = FALSE)
        ## Check if the locally cached version wasn't modified.
        if (!is.null(cache) && cache$info$mtime == info$mtime) {
            message("Using GPL object found in memory from locally cached version.")
            return(cache$gpl)
        }
    }
    txt = data.table::fread(fname, sep = "")[[1]]
    # read_lines reads separate blank lines on windows, so remove them before
    # proceeding. NOT doing so results in the Table header being repeated as the
    # first line of the Table and test failures galore.
    txt = txt[txt != ""]
    return(.parseGPLTxt(txt))
}


txtGrab <- function(regex, x) {
    x <- as.character(x)
    a <- regexpr(regex, x, perl = TRUE)
    return(substr(x, a, a + attr(a, "match.length") - 1))
}

### Function wrapper to get and parse ALL the GSEMatrix files associated with a
### GSE into a list of ExpressionSets
getAndParseGSEMatrices <- function(GEO, destdir, AnnotGPL, getGPL = TRUE, parseCharacteristics = TRUE) {
    GEO <- toupper(GEO)
    ## This stuff functions to get the listing of available files for a given GSE
    ## given that there may be many GSEMatrix files for a given GSE.
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
    b = getDirListing(sprintf(gdsurl, stub, GEO))
    message(sprintf("Found %d file(s)", length(b)))
    ret <- list()
    ## Loop over the files, returning a list, one element for each file
    for (i in 1:length(b)) {
        message(b[i])
        destfile = file.path(destdir, b[i])
        if (file.exists(destfile)) {
            message(sprintf("Using locally cached version: %s", destfile))
        } else {
            url = sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s",
                stub, GEO, b[i])
            downloadFile(url, destfile = destfile, mode = "wb")
        }
        ret[[b[i]]] <- parseGSEMatrix(destfile, destdir = destdir, AnnotGPL = AnnotGPL,
            getGPL = getGPL)$eset
    }
    return(ret)
}

#' Parse a GSE mstrix file
#'
#' Not meant for user calling, but parses
#' a single GSEMatrix file.
#'
#' @importFrom dplyr select filter mutate mutate_all
#' @importFrom tidyr gather spread separate
#' @importFrom readr read_lines
#' @importClassesFrom Biobase ExpressionSet
#' @importFrom magrittr %>%
#'
#' @param fname filename
#' @param AnnotGPL set to TRUE to get the annotation GPL version
#' @param destdir the destdination directory for download
#' @param getGPL whether or not to get the GPL associated
#' @param parseCharacteristics Whether or not to do full 'characteristic' parsing
#' @keywords internal
#' 
parseGSEMatrix <- function(fname, AnnotGPL = FALSE, destdir = tempdir(), getGPL = TRUE,
    parseCharacteristics = TRUE) {
    # Convenient VERY fast line reader based on data.table See:
    # https://stackoverflow.com/a/32924981/459633 data read into single character
    # column, so subset to get just text.
    dat <- data.table::fread(fname, sep = "")[[1]]

    ## get the number of !Series and !Sample lines
    series_header_row_count <- sum(grepl("^!Series_", dat))
    # In the case of ^M in the metadata (GSE781, for example), the line counts
    # for 'skip' and read.table are different.  This next line gets the 'skip'
    # line count for use below in the tmpdat reading 'skip'
    sample_header_start <- grep("^!Sample_", dat)[1]
    samples_header_row_count <- sum(grepl("^!Sample_", dat))
    series_table_begin_line = grep("^!series_matrix_table_begin", dat)
    series_table_end_line = grep("^!series_matrix_table_end", dat)
    if (length(series_table_begin_line) != 1) {
        stop("parsing failed--expected only one '!series_data_table_begin'")
    }
    # con <- fileOpen(fname) Read the !Series_ and !Sample_ lines
    header <- data.table::fread(fname, header = FALSE, nrows = series_header_row_count)
    tmpdat <- data.table::fread(fname, header = FALSE, nrows = samples_header_row_count,
        skip = sample_header_start - 1)

    headertmp <- t(header)
    headerdata <- rbind(data.frame(), headertmp[-1, ])
    colnames(headerdata) <- sub("!Series_", "", as.character(header[[1]]))
    headerlist <- lapply(split.default(headerdata, names(headerdata)), function(x) {
        as.character(Reduce(function(a, b) {
            paste(a, b, sep = "\n")
        }, x))
    })

    link = "https://www.ncbi.nlm.nih.gov/geo/"
    if (!is.null(headerlist$web_link)) {
        link <- headerlist$web_link
    } else if (!is.null(headerlist$geo_accession)) {
        link <- paste(link, "query/acc.cgi?acc=", headerlist$geo_accession, sep = "")
    }

    ed <- new("MIAME", name = ifelse(is.null(headerlist$contact_name), "", headerlist$contact_name),
        title = headerlist$title, contact = ifelse(is.null(headerlist$contact_email),
            "", headerlist$contact_email), pubMedIds = ifelse(is.null(headerlist$pubmed_id),
            "", headerlist$pubmed_id), abstract = ifelse(is.null(headerlist$summary),
            "", headerlist$summary), url = link, other = headerlist)

    tmptmp <- t(tmpdat)
    sampledat <- rbind(data.frame(), tmptmp[-1, ])
    colnames(sampledat) <- make.unique(sub("!Sample_", "", as.character(tmpdat[[1]])))
    sampledat[["geo_accession"]] = as.character(sampledat[["geo_accession"]])
    rownames(sampledat) = sampledat[["geo_accession"]]
    ## Lots of GSEs now use 'characteristics_ch1' and 'characteristics_ch2' for
    ## key-value pairs of annotation. If that is the case, this simply cleans
    ## those up and transforms the keys to column names and the values to column
    ## values.
    if (length(grep("characteristics_ch", colnames(sampledat))) > 0 && parseCharacteristics) {
        pd = sampledat %>%
            dplyr::select(dplyr::contains("characteristics_ch")) %>%
            dplyr::mutate(accession = rownames(.)) %>%
            # these next two lines avoid warnings due to columns having different
            # factor levels (attributes).
        mutate_all(as.character) %>%
            tidyr::gather(characteristics, kvpair, -accession) %>%
            dplyr::filter(grepl(":", kvpair) & !is.na(kvpair))
        # Thx to Mike Smith (@grimbough) for this code sometimes the
        # 'characteristics_ch1' fields are empty and contain no key:value pairs.
        # spread() will fail when called on an empty data_frame.  We catch this
        # case and remove the 'charactics_ch1' column instead
        if (nrow(pd)) {
            pd = dplyr::mutate(pd, characteristics = ifelse(grepl("_ch2", characteristics),
                "ch2", "ch1")) %>%
                tidyr::separate(kvpair, into = c("k", "v"), sep = ":", fill = "right",
                  extra = "merge") %>%
                dplyr::mutate(k = paste(k, characteristics, sep = ":")) %>%
                dplyr::select(-characteristics) %>%
                dplyr::filter(!is.na(v)) %>%
                dplyr::group_by(accession, k) %>%
                dplyr::mutate(v = paste0(trimws(v), collapse = ";")) %>%
                unique() %>%
                tidyr::spread(k, v)
        } else {
            pd = pd %>%
                dplyr::select(accession)
        }
        ## dplyr::mutate(characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1'))
        ## %>% dplyr::filter(grepl(':',kvpair)) %>% tidyr::separate(kvpair, into=
        ## c('k','v'), sep=':') if(nrow(pd)>0) { pd = pd %>% dplyr::mutate(k =
        ## paste(k,characteristics,sep=':')) %>% dplyr::select(-characteristics)
        ## %>% tidyr::spread(k,v)
        sampledat = sampledat %>%
            dplyr::left_join(pd, by = c(geo_accession = "accession"))
    }

    datamat = NULL
    if (series_table_end_line - series_table_begin_line == 2) {
        datamat = read.table(textConnection(dat[(series_table_begin_line + 1):(series_table_end_line -
            1)]), header = TRUE, sep = "\t")
    } else {
        datamat <- data.table::fread(text = dat[(series_table_begin_line + 1):(series_table_end_line -
            1)], quote = "\"", na.strings = c("NA", "null", "NULL", "Null"))
    }
    ## kip = series_table_begin_line) comment.char = '!series_matrix_table_end')
    tmprownames = as.character(datamat[[1]])
    # need the as.matrix for single-sample or empty GSE
    datamat <- as.matrix(datamat[!is.na(tmprownames), -1])
    rownames(datamat) <- tmprownames[!is.na(tmprownames)]
    datamat <- as.matrix(datamat)
    rownames(sampledat) <- colnames(datamat)
    GPL = as.character(sampledat[1, grep("^platform_id", colnames(sampledat), ignore.case = TRUE)])
    ## if getGPL is FALSE, skip this and featureData is then a data.frame with no
    ## columns
    fd = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(datamat)))
    if (getGPL) {
        gpl <- getGEO(GPL, AnnotGPL = AnnotGPL, destdir = destdir)
        vmd <- Columns(gpl)
        dat <- Table(gpl)
        ## GEO uses 'TAG' instead of 'ID' for SAGE GSE/GPL entries, but it is
        ## apparently always the first column, so use dat[,1] instead of dat$ID
        ## The next line deals with the empty GSE
        tmpnames = character(0)
        if (ncol(dat) > 0) {
            tmpnames = as.character(dat[, 1])
        }
        ## Fixed bug caused by an ID being 'NA' in GSE15197, for example
        tmpnames[is.na(tmpnames)] = "NA"
        rownames(dat) <- make.unique(tmpnames)
        ## Apparently, NCBI GEO uses case-insensitive matching between platform
        ## IDs and series ID Refs ???
        dat <- dat[match(tolower(rownames(datamat)), tolower(rownames(dat))), ]
        # Fix possibility of duplicate column names in the GPL files; this is
        # prevalent in the Annotation GPLs
        rownames(vmd) <- make.unique(colnames(Table(gpl)))
        colnames(dat) <- rownames(vmd)
        fd <- new("AnnotatedDataFrame", data = dat, varMetadata = vmd)
    }
    if (is.null(nrow(datamat))) {
        ## fix empty GSE datamatrix samplename stuff above does not work with
        ## empty GSEs, so fix here, also
        tmpnames <- names(datamat)
        rownames(sampledat) <- tmpnames
        datamat = matrix(nrow = 0, ncol = nrow(sampledat))
        colnames(datamat) <- tmpnames
    } else {
        ## This looks like a dangerous operation but is needed to deal with the
        ## fact that NCBI GEO allows case-insensitive matching and we need to
        ## pick one.
        rownames(datamat) <- rownames(dat)
    }
    eset <- new("ExpressionSet", phenoData = as(sampledat, "AnnotatedDataFrame"), annotation = GPL,
        featureData = fd, experimentData = ed, exprs = as.matrix(datamat))
    return(list(GPL = as.character(sampledat[1, grep("platform_id", colnames(sampledat),
        ignore.case = TRUE)]), eset = eset))
}
