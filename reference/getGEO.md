# Get a GEO object from NCBI or file

This function is the main user-level function in the GEOquery package.
It directs the download (if no filename is specified) and parsing of a
GEO SOFT format file into an R data structure specifically designed to
make access to each of the important parts of the GEO SOFT format easily
accessible.

## Usage

``` r
getGEO(
  GEO = NULL,
  filename = NULL,
  destdir = tempdir(),
  GSElimits = NULL,
  GSEMatrix = TRUE,
  AnnotGPL = FALSE,
  getGPL = TRUE,
  parseCharacteristics = TRUE
)
```

## Arguments

- GEO:

  A character string representing a GEO object for download and parsing.
  (eg., 'GDS505','GSE2','GSM2','GPL96')

- filename:

  The filename of a previously downloaded GEO SOFT format file or its
  gzipped representation (in which case the filename must end in .gz).
  Either one of GEO or filename may be specified, not both. GEO series
  matrix files are also handled. Note that since a single file is being
  parsed, the return value is not a list of esets, but a single eset
  when GSE matrix files are parsed.

- destdir:

  The destination directory for any downloads. Defaults to the
  architecture-dependent tempdir. You may want to specify a different
  directory if you want to save the file for later use. Doing so is a
  good idea if you have a slow connection, as some of the GEO files are
  HUGE!

- GSElimits:

  This argument can be used to load only a contiguous subset of the GSMs
  from a GSE. It should be specified as a vector of length 2 specifying
  the start and end (inclusive) GSMs to load. This could be useful for
  splitting up large GSEs into more manageable parts, for example.

- GSEMatrix:

  A boolean telling GEOquery whether or not to use GSE Series Matrix
  files from GEO. The parsing of these files can be many
  orders-of-magnitude faster than parsing the GSE SOFT format files.
  Defaults to TRUE, meaning that the SOFT format parsing will not occur;
  set to FALSE if you for some reason need other columns from the GSE
  records.

- AnnotGPL:

  A boolean defaulting to FALSE as to whether or not to use the
  Annotation GPL information. These files are nice to use because they
  contain up-to-date information remapped from Entrez Gene on a regular
  basis. However, they do not exist for all GPLs; in general, they are
  only available for GPLs referenced by a GDS

- getGPL:

  A boolean defaulting to TRUE as to whether or not to download and
  include GPL information when getting a GSEMatrix file. You may want to
  set this to FALSE if you know that you are going to annotate your
  featureData using Bioconductor tools rather than relying on
  information provided through NCBI GEO. Download times can also be
  greatly reduced by specifying FALSE.

- parseCharacteristics:

  A boolean defaulting to TRUE as to whether or not to parse the
  characteristics information (if available) for a GSE Matrix file. Set
  this to FALSE if you experience trouble while parsing the
  characteristics.

## Value

An object of the appropriate class (GDS, GPL, GSM, or GSE) is returned.
If the GSEMatrix option is used, then a list of ExpressionSet objects is
returned, one for each SeriesMatrix file associated with the GSE
accesion. If the filename argument is used in combination with a
GSEMatrix file, then the return value is a single ExpressionSet.

## Details

getGEO functions to download and parse information available from NCBI
GEO (<http://www.ncbi.nlm.nih.gov/geo>). Here are some details about
what is avaible from GEO. All entity types are handled by getGEO and
essentially any information in the GEO SOFT format is reflected in the
resulting data structure.

From the GEO website:

The Gene Expression Omnibus (GEO) from NCBI serves as a public
repository for a wide range of high-throughput experimental data. These
data include single and dual channel microarray-based experiments
measuring mRNA, genomic DNA, and protein abundance, as well as non-array
techniques such as serial analysis of gene expression (SAGE), and mass
spectrometry proteomic data. At the most basic level of organization of
GEO, there are three entity types that may be supplied by users:
Platforms, Samples, and Series. Additionally, there is a curated entity
called a GEO dataset.

A Platform record describes the list of elements on the array (e.g.,
cDNAs, oligonucleotide probesets, ORFs, antibodies) or the list of
elements that may be detected and quantified in that experiment (e.g.,
SAGE tags, peptides). Each Platform record is assigned a unique and
stable GEO accession number (GPLxxx). A Platform may reference many
Samples that have been submitted by multiple submitters.

A Sample record describes the conditions under which an individual
Sample was handled, the manipulations it underwent, and the abundance
measurement of each element derived from it. Each Sample record is
assigned a unique and stable GEO accession number (GSMxxx). A Sample
entity must reference only one Platform and may be included in multiple
Series.

A Series record defines a set of related Samples considered to be part
of a group, how the Samples are related, and if and how they are
ordered. A Series provides a focal point and description of the
experiment as a whole. Series records may also contain tables describing
extracted data, summary conclusions, or analyses. Each Series record is
assigned a unique and stable GEO accession number (GSExxx).

GEO DataSets (GDSxxx) are curated sets of GEO Sample data. A GDS record
represents a collection of biologically and statistically comparable GEO
Samples and forms the basis of GEO's suite of data display and analysis
tools. Samples within a GDS refer to the same Platform, that is, they
share a common set of probe elements. Value measurements for each Sample
within a GDS are assumed to be calculated in an equivalent manner, that
is, considerations such as background processing and normalization are
consistent across the dataset. Information reflecting experimental
design is provided through GDS subsets.

## Warning

Some of the files that are downloaded, particularly those associated
with GSE entries from GEO are absolutely ENORMOUS and parsing them can
take quite some time and memory. So, particularly when working with
large GSE entries, expect that you may need a good chunk of memory and
that coffee may be involved when parsing....

## See also

[`getGEOfile`](http://seandavi.github.io/GEOquery/reference/getGEOfile.md)

## Author

Sean Davis

## Examples

``` r

gds <- getGEO('GDS10')
gds
#> An object of class "GDS"
#> channel_count 
#> [1] "1"
#> dataset_id 
#>  [1] "GDS10" "GDS10" "GDS10" "GDS10" "GDS10" "GDS10" "GDS10" "GDS10" "GDS10"
#> [10] "GDS10" "GDS10" "GDS10"
#> description 
#>  [1] "Examination of spleen and thymus of type 1 diabetes nonobese diabetic (NOD) mouse, four NOD-derived diabetes-resistant congenic strains and two nondiabetic control strains."
#>  [2] "spleen"                                                                                                                                                                      
#>  [3] "thymus"                                                                                                                                                                      
#>  [4] "NOD"                                                                                                                                                                         
#>  [5] "Idd3"                                                                                                                                                                        
#>  [6] "Idd5"                                                                                                                                                                        
#>  [7] "Idd3+Idd5"                                                                                                                                                                   
#>  [8] "Idd9"                                                                                                                                                                        
#>  [9] "B10.H2g7"                                                                                                                                                                    
#> [10] "B10.H2g7 Idd3"                                                                                                                                                               
#> [11] "diabetic"                                                                                                                                                                    
#> [12] "diabetic-resistant"                                                                                                                                                          
#> [13] "nondiabetic"                                                                                                                                                                 
#> email 
#> [1] "geo@ncbi.nlm.nih.gov"
#> feature_count 
#> [1] "39114"
#> institute 
#> [1] "NCBI NLM NIH"
#> name 
#> [1] "Gene Expression Omnibus (GEO)"
#> order 
#> [1] "none"
#> platform 
#> [1] "GPL24"
#> platform_organism 
#> [1] "Mus musculus"
#> platform_technology_type 
#> [1] "in situ oligonucleotide"
#> pubmed_id 
#> [1] "11827943"
#> ref 
#> [1] "Nucleic Acids Res. 2005 Jan 1;33 Database Issue:D562-6"
#> reference_series 
#> [1] "GSE11"
#> sample_count 
#> [1] "28"
#> sample_id 
#>  [1] "GSM582,GSM583,GSM584,GSM585,GSM586,GSM587,GSM588,GSM589,GSM590,GSM591,GSM592,GSM593,GSM594,GSM595"              
#>  [2] "GSM596,GSM597,GSM598,GSM599,GSM600,GSM601,GSM602,GSM603,GSM604,GSM605,GSM606,GSM607,GSM608,GSM609"              
#>  [3] "GSM582,GSM589,GSM596,GSM603"                                                                                    
#>  [4] "GSM583,GSM590,GSM597,GSM604"                                                                                    
#>  [5] "GSM584,GSM591,GSM598,GSM605"                                                                                    
#>  [6] "GSM585,GSM592,GSM599,GSM606"                                                                                    
#>  [7] "GSM586,GSM593,GSM600,GSM607"                                                                                    
#>  [8] "GSM587,GSM594,GSM601,GSM608"                                                                                    
#>  [9] "GSM588,GSM595,GSM602,GSM609"                                                                                    
#> [10] "GSM582,GSM589,GSM596,GSM603"                                                                                    
#> [11] "GSM583,GSM590,GSM597,GSM604,GSM584,GSM591,GSM598,GSM605,GSM585,GSM592,GSM599,GSM606,GSM586,GSM593,GSM600,GSM607"
#> [12] "GSM587,GSM594,GSM601,GSM608,GSM588,GSM595,GSM602,GSM609"                                                        
#> sample_organism 
#> [1] "Mus musculus"
#> sample_type 
#> [1] "RNA"
#> title 
#> [1] "Type 1 diabetes gene expression profiling"
#> type 
#>  [1] "Expression profiling by array" "tissue"                       
#>  [3] "tissue"                        "strain"                       
#>  [5] "strain"                        "strain"                       
#>  [7] "strain"                        "strain"                       
#>  [9] "strain"                        "strain"                       
#> [11] "disease state"                 "disease state"                
#> [13] "disease state"                
#> update_date 
#> [1] "Jul 15 2003"
#> value_type 
#> [1] "count"
#> web_link 
#> [1] "http://www.ncbi.nlm.nih.gov/geo"
#> An object of class "GEODataTable"
#> ****** Column Descriptions ******
#>    sample tissue        strain      disease.state
#> 1  GSM582 spleen           NOD           diabetic
#> 2  GSM589 spleen           NOD           diabetic
#> 3  GSM583 spleen          Idd3 diabetic-resistant
#> 4  GSM590 spleen          Idd3 diabetic-resistant
#> 5  GSM584 spleen          Idd5 diabetic-resistant
#> 6  GSM591 spleen          Idd5 diabetic-resistant
#> 7  GSM585 spleen     Idd3+Idd5 diabetic-resistant
#> 8  GSM592 spleen     Idd3+Idd5 diabetic-resistant
#> 9  GSM586 spleen          Idd9 diabetic-resistant
#> 10 GSM593 spleen          Idd9 diabetic-resistant
#> 11 GSM587 spleen      B10.H2g7        nondiabetic
#> 12 GSM594 spleen      B10.H2g7        nondiabetic
#> 13 GSM588 spleen B10.H2g7 Idd3        nondiabetic
#> 14 GSM595 spleen B10.H2g7 Idd3        nondiabetic
#> 15 GSM596 thymus           NOD           diabetic
#> 16 GSM603 thymus           NOD           diabetic
#> 17 GSM597 thymus          Idd3 diabetic-resistant
#> 18 GSM604 thymus          Idd3 diabetic-resistant
#> 19 GSM598 thymus          Idd5 diabetic-resistant
#> 20 GSM605 thymus          Idd5 diabetic-resistant
#> 21 GSM599 thymus     Idd3+Idd5 diabetic-resistant
#> 22 GSM606 thymus     Idd3+Idd5 diabetic-resistant
#> 23 GSM600 thymus          Idd9 diabetic-resistant
#> 24 GSM607 thymus          Idd9 diabetic-resistant
#> 25 GSM601 thymus      B10.H2g7        nondiabetic
#> 26 GSM608 thymus      B10.H2g7        nondiabetic
#> 27 GSM602 thymus B10.H2g7 Idd3        nondiabetic
#> 28 GSM609 thymus B10.H2g7 Idd3        nondiabetic
#>                                        description
#> 1            Value for GSM582: NOD_S1; src: Spleen
#> 2            Value for GSM589: NOD_S2; src: Spleen
#> 3           Value for GSM583: Idd3_S1; src: Spleen
#> 4           Value for GSM590: Idd3_S2; src: Spleen
#> 5           Value for GSM584: Idd5_S1; src: Spleen
#> 6           Value for GSM591: Idd5_S2; src: Spleen
#> 7         Value for GSM585: Idd3+5_S1; src: Spleen
#> 8         Value for GSM592: Idd3+5_S2; src: Spleen
#> 9           Value for GSM586: Idd9_S1; src: Spleen
#> 10          Value for GSM593: Idd9_S2; src: Spleen
#> 11      Value for GSM587: B10.H2g7_S1; src: Spleen
#> 12      Value for GSM594: B10.H2g7_S2; src: Spleen
#> 13 Value for GSM588: B10.H2g7 Idd3_S1; src: Spleen
#> 14 Value for GSM595: B10.H2g7 Idd3_S2; src: Spleen
#> 15           Value for GSM596: NOD_T1; src: Thymus
#> 16           Value for GSM603: NOD_T2; src: Thymus
#> 17          Value for GSM597: Idd3_T1; src: Thymus
#> 18          Value for GSM604: Idd3_T2; src: Thymus
#> 19          Value for GSM598: Idd5_T1; src: Thymus
#> 20          Value for GSM605: Idd5_T2; src: Thymus
#> 21        Value for GSM599: Idd3+5_T1; src: Thymus
#> 22        Value for GSM606: Idd3+5_T2; src: Thymus
#> 23          Value for GSM600: Idd9_T1; src: Thymus
#> 24          Value for GSM607: Idd9_T2; src: Thymus
#> 25      Value for GSM601: B10.H2g7_T1; src: Thymus
#> 26      Value for GSM608: B10.H2g7_T2; src: Thymus
#> 27 Value for GSM602: B10.H2g7 Idd3_T1; src: Thymus
#> 28 Value for GSM609: B10.H2g7 Idd3_T2; src: Thymus
#> ****** Data Table ******

gse <- getGEO('GSE10')
#> Found 1 file(s)
#> GSE10_series_matrix.txt.gz
# Returns a list, so look at first item

gse[[1]]
#> ExpressionSet (storageMode: lockedEnvironment)
#> assayData: 96903 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: GSM571 GSM572 GSM573 GSM574
#>   varLabels: title geo_accession ... data_row_count (34 total)
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: AAAAAAAAAA AAAAAAAAAC ... TTTTTTTTTT (96903 total)
#>   fvarLabels: TAG GI
#>   fvarMetadata: Column Description labelDescription
#> experimentData: use 'experimentData(object)'
#>   pubMedIds: 11756676 
#> Annotation: GPL4 
```
