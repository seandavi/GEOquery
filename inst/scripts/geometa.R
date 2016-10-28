require(rentrez)
require(purrr)

.docSumListConvert = function(x) {
    ret = data.frame(
        uid          = x %>% purrr::map('uid') %>% purrr::flatten_chr() %>% as.integer(),
        accession    = x %>% purrr::map('accession') %>% purrr::flatten_chr(),
        seriestitle  = x %>% purrr::map('seriestitle') %>% purrr::flatten_chr(),
        platformtitle= x %>% purrr::map('platformtitle') %>% purrr::flatten_chr(),
        taxon        = x %>% purrr::map('gpl') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% I(),
        entrytype    = x %>% purrr::map('entrytype') %>% purrr::flatten_chr(),
        gdstype      = x %>% purrr::map('gdstype') %>% purrr::flatten_chr(),
        ptechtype    = x %>% purrr::map('ptechtype') %>% purrr::flatten_chr(),
        valtype      = x %>% purrr::map('valtype') %>% purrr::flatten_chr(),
        ssinfo       = x %>% purrr::map('ssinfo') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% I(),
        title        = x %>% purrr::map('title') %>% purrr::flatten_chr(),
        summary      = x %>% purrr::map('summary') %>% purrr::flatten_chr(),
        gpl          = x %>% purrr::map('gpl') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% purrr::map(function(x) paste0('GPL',x)) %>% I(),
        gse          = x %>% purrr::map('gse') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% purrr::map(function(x) paste0('GSE',x)) %>% I(),
        gds          = x %>% purrr::map('gds') %>% purrr::map(strsplit,';') %>% purrr::flatten() %>% purrr::map(function(x) ifelse(nchar(x)>0,paste0('GDS',x),"")) %>% I(),
        samples      = x %>% purrr::map('samples') %>% purrr::map('accession') %>% I(),
        suppfile     = x %>% purrr::map('suppfile') %>% purrr::map(strsplit,', ') %>% purrr::flatten() %>% I(),
        ftplink      = x %>% purrr::map('ftplink') %>% purrr::flatten_chr(),
        n_samples    = x %>% purrr::map('n_samples') %>% as.integer(),
        pubmedids    = x %>% purrr::map('pubmedids') %>% I(),
        projects     = x %>% purrr::map('projects') %>% I(),
        public_date  = x %>% purrr::map('pdat') %>% purrr::flatten_chr() %>% date(),
        samplestaxa  = x %>% purrr::map('samplestaxa') %>% I() #purrr::map(strsplit,', ') %>% purrr::flatten() %>% I()
        )
    return(ret)
    }

#' get all GEO records as docsums
#'
#' Fetches all records from NCBI entrez and creates a data frame
#' from the docsum entries, loadable into a SQL-like database
#'
#'
searchGEODocSums = function(term='GPL[ETYP] OR GSE[ETYP] OR GSM[ETYP] OR GDS[ETYP]') {
    res = rentrez::entrez_search(db='gds',use_history=TRUE,term=term)
    return(res)
}
fetchDocSums= function(res,retstart=1,retmax=100) {
  z=NULL
  while(is.null(z) || inherits(z,'error')) {
    z = try(rentrez::entrez_summary(db='gds',
                                    web_history=res$web_history,
                                    retstart=retstart,
                                    retmax=retmax))
  }
    
  return(.docSumListConvert(z))
}




getGEOMeta = function(geo) {
    .getGSMmeta = function(x) {
        ret = list(
            accession     = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Sample'),'iid'),
            gpl           = xml_attr(xml_find_first(dat,'/d1:MINiML/d1:Platform'),'iid'),
            title         = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Title')),
            type          = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Type')),
            submission_date= xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Status/d1:Submission-Date')),
            last_update   = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Status/d1:Last-Update--Date')),
            release_date  = xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Status/d1:Release-Date')),
            n_channels    = as.integer(xml_text(xml_find_first(dat,'/d1:MINiML/d1:Sample/d1:Channel-Count'))),
            channels      = data.frame(do.call(rbind,lapply(xml_find_all(dat,'/d1:MINiML/d1:Sample/d1:Channel'),
                                   function(channel) {
                                       return(list(
                                           tax_id = xml_attr(xml_find_first(channel,'./d1:Organism'),'taxid'),
                                           organixm = xml_text(xml_find_first(channel,'./d1:Organism')),
                                           source = xml_text(xml_find_first(channel,'./d1:Source')),
                                           molecule = xml_text(xml_find_first(channel,'./d1:Molecule')),
                                           characteristics = xml_text(xml_find_first(channel,'./d1:Characteristics'))
                                           
                                       ))
                                   })
                               
        )))
        return(ret)
    }
    geo = toupper(geo)
    stopifnot(substr(geo,1,3) %in% c('GSE','GSM','GPL','GDS'))
              
    dat = read_xml(sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=full&targ=self&form=xml",geo))
    ret = switch(substr(geo,1,3),
                 GSM = .getGSMmeta(dat),
                 NA)
    return(ret)
}    
