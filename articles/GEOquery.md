# Using the GEOquery Package

## Introduction to GEO and GEOquery

The NCBI Gene Expression Omnibus
([GEO](https://www.ncbi.nlm.nih.gov/geo/)) was established in 2000 as a
public repository for high-throughput molecular abundance data,
primarily microarray data at the time. Today, GEO hosts a diverse array
of data types including gene expression, genomic DNA, and protein
abundance measurements from various technologies like microarrays,
next-generation sequencing, and mass spectrometry.

GEOquery was created to bridge the gap between this vast resource and
Bioconductor’s analytical tools. First released in 2007, GEOquery has
evolved alongside GEO itself, adapting to new data types and formats
over time.

### Why GEOquery?

Before GEOquery, researchers would need to:

1.  Manually download data from the GEO website
2.  Parse complex SOFT format files
3.  Construct data structures suitable for analysis
4.  Integrate metadata with expression data

GEOquery automates this entire process, allowing researchers to focus on
analysis rather than data acquisition and formatting. GEOquery also
facilitates automation and reproducibility by incorporating data
acquisition into workflows, scripts, or documents.

### GEO Data Organization

Understanding GEO’s data organization is essential for effective use of
GEOquery:

1.  **Platform (GPL)**: Describes array design, probes, or detectable
    elements
2.  **Sample (GSM)**: Contains individual experiment measurements
3.  **Series (GSE)**: Groups related samples together, typically
    representing a complete study
4.  **Dataset (GDS)**: Curated by GEO staff, represents biologically and
    statistically comparable samples

## Getting Started with GEOquery

Before working with GEOquery (or any R or Bioconductor package), one
must first install the package. Installation is a one-time (or at least
not repeated often) operation. To install GEOquery, ensure that you have
[installed R and Bioconductor](https://bioconductor.org/install/).

Then, to install GEOquery:

``` r
BiocManager::install('GEOquery')
```

Before using GEOquery, we need to load the GEOquery library. Loading the
GEOquery library must be done *each time you start a new R session*.

``` r
library(GEOquery)
```

### Downloading a GEO Series

The most common use case is downloading a GEO Series (GSE), which
typically represents a complete study:

``` r
# Download GSE2553
gse <- getGEO("GSE2553")
```

    Found 1 file(s)

    GSE2553_series_matrix.txt.gz

``` r
class(gse)
```

    [1] "list"

Notice that `getGEO` returns a list. This is because a single GSE can
contain experiments from multiple platforms. Each element of the list is
an `ExpressionSet` containing data from one platform:

``` r
length(gse)
```

    [1] 1

``` r
gse[[1]]
```

    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 12600 features, 181 samples
      element names: exprs
    protocolData: none
    phenoData
      sampleNames: GSM48681 GSM48682 ... GSM48861 (181 total)
      varLabels: title geo_accession ... data_row_count (30 total)
      varMetadata: labelDescription
    featureData
      featureNames: 1 2 ... 12600 (12600 total)
      fvarLabels: ID PenAt ... Chimeric_Cluster_IDs (13 total)
      fvarMetadata: Column Description labelDescription
    experimentData: use 'experimentData(object)'
      pubMedIds: 16230383
    Annotation: GPL1977 

### Historical Context: SOFT Format vs GSEMatrix Files

GEO originally provided data in SOFT (Simple Omnibus Format in Text)
format, which contained extensive information but was slow to parse for
large datasets. In response to community needs, GEO introduced GSEMatrix
files—a more efficient, tab-delimited format.

GEOquery defaults to using GSEMatrix files (`GSEMatrix=TRUE`)
because: 1. Parsing is substantially faster (often by 10-100x) 2. Memory
usage is more efficient 3. The resulting `ExpressionSet` objects are
directly usable with Bioconductor tools

## Searching GEO Programmatically

While GEO’s web interface is powerful, programmatic searches enable
automated data discovery and retrieval. GEOquery provides direct access
to GEO’s search capabilities:

``` r
# What fields can we search?
fields <- searchFieldsGEO()
kable(fields)
```

| Name | FullName     | Description | TermCount | IsDate | IsNumerical | SingleToken | Hierarchy | IsHidden |
|:-----|:-------------|:------------|:----------|:-------|:------------|:------------|:----------|:---------|
| ALL  | All Fields   | All term….  | 45152635  | N      | N           | N           | N         | N        |
| UID  | UID          | Unique n….  | 0         | N      | Y           | Y           | N         | Y        |
| FILT | Filter       | Limits t….  | 71        | N      | N           | Y           | N         | N        |
| ORGN | Organism     | exploded….  | 75623     | N      | N           | Y           | Y         | N        |
| ACCN | GEO Acce….   | accessio….  | 19303976  | N      | N           | Y           | N         | N        |
| TITL | Title        | Words in….  | 9946487   | N      | N           | Y           | N         | N        |
| DESC | Description  | Text fro….  | 10548780  | N      | N           | Y           | N         | N        |
| SFIL | Suppleme….   | Suppleme….  | 255       | N      | N           | Y           | N         | N        |
| ETYP | Entry Type   | Entry ty….  | 4         | N      | N           | Y           | N         | N        |
| STYP | Sample Type  | Sample type | 9         | N      | N           | Y           | N         | N        |
| VTYP | Sample V….   | type of ….  | 7         | N      | N           | Y           | N         | N        |
| PTYP | Platform….   | Platform….  | 17        | N      | N           | Y           | N         | N        |
| GTYP | DataSet Type | type of ….  | 27        | N      | N           | Y           | N         | N        |
| NSAM | Number o….   | Number o….  | 2134      | N      | Y           | Y           | N         | N        |
| SRC  | Sample S….   | sample s….  | 461250    | N      | N           | Y           | N         | N        |
| AUTH | Author       | author o….  | 1300542   | N      | N           | Y           | N         | N        |
| INST | Submitte….   | institut….  | 25314     | N      | N           | Y           | N         | N        |
| NPRO | Number o….   | number o….  | 7258      | N      | Y           | Y           | N         | N        |
| SSTP | Subset V….   | subset v….  | 24        | N      | N           | Y           | N         | N        |
| SSDE | Subset D….   | subset d….  | 7535      | N      | N           | Y           | N         | N        |
| GEID | Reporter….   | name or ….  | 2840498   | N      | N           | Y           | N         | N        |
| PDAT | Publicat….   | publicat….  | 8453      | Y      | N           | Y           | N         | N        |
| UDAT | Update Date  | date        | 7570      | Y      | N           | Y           | N         | N        |
| TAGL | Tag Length   | Tag/Sign….  | 9         | N      | N           | Y           | N         | N        |
| RGSE | Related ….   | Related ….  | 30135     | N      | N           | Y           | N         | N        |
| RGPL | Related ….   | Related ….  | 271071    | N      | N           | Y           | N         | N        |
| MESH | MeSH Terms   | Medical ….  | 17643     | N      | N           | Y           | Y         | N        |
| PROJ | Project      | Project     | 10        | N      | N           | Y           | N         | N        |
| ATNM | Attribut….   | Attribut….  | 49222     | N      | N           | Y           | N         | N        |
| ATTR | Attribute    | Attribute   | 2763893   | N      | N           | Y           | N         | N        |
| PROP | Properties   | Properties  | 3         | N      | N           | Y           | N         | N        |

Table 1: Available GEO search fields.

GEO uses a specific search syntax with field identifiers in square
brackets:

``` r
# Find RNA-seq studies related to COVID-19 in humans
results <- searchGEO('covid-19[All Fields] AND "rnaseq counts"[Filter] AND Homo sapiens[ORGN]')
results |>
  dplyr::mutate(Summary=paste(strtrim(Summary,120), '...')) |>
  dplyr::mutate(Title = paste(strtrim(Title, 120), '...')) |>
  head() |>
  kable()
```

| Title                                                                                                                      | Summary                                                                                                                    | Organism     | Type                                               | Platforms | Contains   | FTP download                                                              | Series Accession |        ID | SRA Run Selector |
|:---------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------|:-------------|:---------------------------------------------------|:----------|:-----------|:--------------------------------------------------------------------------|:-----------------|----------:|:-----------------|
| Comparative Analysis of Host Responses Between Sepsis and COVID-19: A Prospective Observational Study of Whole Blood Tra … | This prospective observational study conducted at Osaka University Graduate School of Medicine aimed to compare host res … | Homo sapiens | Expression profiling by high throughput sequencing | GPL30209  | 72 Samples | GEO (TXT) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243217/      | GSE243217        | 200243217 | NA               |
| An aberrant immune-epithelial progenitor niche drives post-viral lung sequelae \[human\] …                                 | Respiratory viral infections are being increasingly recognized not just for their acute impact but also as potential tri … | Homo sapiens | Other                                              | GPL24676  | 5 Samples  | GEO (H5, PNG) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267226/  | GSE267226        | 200267226 | NA               |
| Longitudinal transcriptomic analysis reveals persistent enrichment of iron homeostasis and erythrocyte function pathways … | The acute respiratory distress syndrome (ARDS) is a common complications of severe COVID-19 and contributes to patient m … | Homo sapiens | Expression profiling by high throughput sequencing | GPL34284  | 49 Samples | GEO (TXT) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273149/      | GSE273149        | 200273149 | NA               |
| Features of chronic urticaria after COVID-19 mRNA vaccine, a real-life cohort study …                                      | New onsets of chronic urticaria (CU) have been reported after repeated immunizations, mainly with the Moderna mRNA-1273 …  | Homo sapiens | Expression profiling by high throughput sequencing | GPL20301  | 32 Samples | GEO (TXT) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272645/      | GSE272645        | 200272645 | NA               |
| Transcriptomic Profiling of Neutrophils and Low-Density Granulocytes in COVID-19 Patients …                                | The severity of COVID-19 is linked to excessive inflammation. Neutrophils represent a critical arm of the innate immune …  | Homo sapiens | Expression profiling by high throughput sequencing | GPL24676  | 36 Samples | GEO (TSV, TXT) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272381/ | GSE272381        | 200272381 | NA               |
| Effects of envelope- or membrane-protein segments of SARS-CoV-2 on gene expression of HUVEC cells …                        | Exterior segments of E-proteins bound onto and modulated gene expression in human vascular endothelial cells in vitro Th … | Homo sapiens | Expression profiling by high throughput sequencing | GPL28038  | 9 Samples  | GEO (TXT) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268369/      | GSE268369        | 200268369 | NA               |

Table 2: Top search results for studies related to COVID-19 in humans
with GEO-calculated RNA-seq counts available.

The search capabilities mirror GEO’s web interface but allow for
integration with R workflows.

## RNA-seq Quantifications: NCBI’s Solution to Reanalysis Challenges

### The Challenge of RNA-seq Reanalysis

A major barrier to exploiting the massive volume of public RNA-seq data
has been the computational cost and expertise required to consistently
process raw reads into usable expression values. Different processing
pipelines can produce different results, making cross-study comparisons
challenging.

### NCBI’s RNA-seq Quantification Pipeline

To address this challenge, in 2020-2021, the NCBI SRA and GEO teams
developed a standardized pipeline that precomputes RNA-seq gene
expression counts for human and mouse datasets. As described in their
[documentation](https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html),
this pipeline:

1.  Processes RNA-seq data from SRA using the HISAT2 aligner
2.  Generates gene expression counts using the featureCounts program
3.  Provides consistent annotation based on current genome builds
4.  Makes counts available in standardized formats

GEOquery provides direct access to these precomputed counts:

``` r
# Check if RNA-seq quantifications are available
has_quant <- hasRNASeqQuantifications("GSE164073")
has_quant
```

    [1] TRUE

``` r
# Get genome build and species information
genome_info <- getRNASeqQuantGenomeInfo("GSE164073")
genome_info
```

``` r
# Download and construct a SummarizedExperiment
se <- getRNASeqData("GSE164073")
se
```

This feature saves researchers significant time and computational
resources while ensuring standardized processing across datasets.

## Understanding Supplementary Files in GEO

GEO accessions often include supplementary files containing raw data,
processing scripts, or additional results not captured in the standard
GEO formats. These files are invaluable for:

1.  Accessing raw data (e.g., CEL files, FASTQ files)
2.  Understanding custom processing pipelines
3.  Retrieving additional metadata or results

GEOquery makes accessing these files straightforward:

``` r
# List available supplementary files without downloading
supp_files <- getGEOSuppFiles('GSE63137', fetch_files = FALSE)
head(supp_files)
```

                                                                 fname
    1                  GSE63137_ATAC-seq_PV_neurons_HOMER_peaks.bed.gz
    2                 GSE63137_ATAC-seq_VIP_neurons_HOMER_peaks.bed.gz
    3          GSE63137_ATAC-seq_excitatory_neurons_HOMER_peaks.bed.gz
    4  GSE63137_ChIP-seq_H3K27ac_excitatory_neurons_SICER_peaks.bed.gz
    5 GSE63137_ChIP-seq_H3K27me3_excitatory_neurons_SICER_peaks.bed.gz
    6  GSE63137_ChIP-seq_H3K4me1_excitatory_neurons_SICER_peaks.bed.gz
                                                                                                                                    url
    1                  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_ATAC-seq_PV_neurons_HOMER_peaks.bed.gz
    2                 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_ATAC-seq_VIP_neurons_HOMER_peaks.bed.gz
    3          https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_ATAC-seq_excitatory_neurons_HOMER_peaks.bed.gz
    4  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_ChIP-seq_H3K27ac_excitatory_neurons_SICER_peaks.bed.gz
    5 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_ChIP-seq_H3K27me3_excitatory_neurons_SICER_peaks.bed.gz
    6  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_ChIP-seq_H3K4me1_excitatory_neurons_SICER_peaks.bed.gz

You can filter files by pattern to find specific file types:

``` r
# Find all text files
txt_files <- getGEOSuppFiles('GSE63137', fetch_files = FALSE, 
                             filter_regex = 'txt')
head(txt_files)
```

                                                         fname
    1                GSE63137_MethylC-seq_DMRs_methylpy.txt.gz
    2         GSE63137_MethylC-seq_PV_neurons_UMRs_LMRs.txt.gz
    3        GSE63137_MethylC-seq_VIP_neurons_UMRs_LMRs.txt.gz
    4 GSE63137_MethylC-seq_excitatory_neurons_UMRs_LMRs.txt.gz
                                                                                                                            url
    1                https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_MethylC-seq_DMRs_methylpy.txt.gz
    2         https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_MethylC-seq_PV_neurons_UMRs_LMRs.txt.gz
    3        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_MethylC-seq_VIP_neurons_UMRs_LMRs.txt.gz
    4 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63137/suppl//GSE63137_MethylC-seq_excitatory_neurons_UMRs_LMRs.txt.gz

And download specific files or all supplementary files:

``` r
# Download all supplementary files for a sample
getGEOSuppFiles('GSM15789') # Files saved to a new directory
```

## Navigating Between R and the GEO Web Interface

Sometimes you may want to examine a GEO record in its web interface.
GEOquery provides convenience functions for this:

``` r
# Get the URL for a GEO accession
url <- urlForAccession("GSE262484")
url

# Open a browser to the GEO page
browseGEOAccession("GSE262484")
```

For RNA-seq datasets specifically, there’s a convenience function to
search for RNA-seq counts on the GEO website:

``` r
browseWebsiteRNASeqSearch()
```

These functions bridge the programmatic and web interfaces to GEO,
allowing seamless transitions between analytical and exploratory modes.

## Working with GDS Datasets

GEO DataSets (GDS) are curated collections of samples, processed and
normalized to be directly comparable. While less common in modern
workflows, they remain available and GEOquery supports them:

``` r
# Download a GDS dataset
gds <- getGEO("GDS507")
gds
```

    An object of class "GDS"
    channel_count
    [1] "1"
    dataset_id
     [1] "GDS507" "GDS507" "GDS507" "GDS507" "GDS507" "GDS507" "GDS507" "GDS507"
     [9] "GDS507" "GDS507" "GDS507" "GDS507"
    description
     [1] "Investigation into mechanisms of renal clear cell carcinogenesis (RCC). Comparison of renal clear cell tumor tissue and adjacent normal tissue isolated from the same surgical samples."
     [2] "RCC"
     [3] "normal"
     [4] "035"
     [5] "023"
     [6] "001"
     [7] "005"
     [8] "011"
     [9] "032"
    [10] "1"
    [11] "2"
    [12] "3"
    [13] "4"
    email
    [1] "geo@ncbi.nlm.nih.gov"
    feature_count
    [1] "22645"
    institute
    [1] "NCBI NLM NIH"
    name
    [1] "Gene Expression Omnibus (GEO)"
    order
    [1] "none"
    platform
    [1] "GPL97"
    platform_organism
    [1] "Homo sapiens"
    platform_technology_type
    [1] "in situ oligonucleotide"
    pubmed_id
    [1] "14641932"
    ref
    [1] "Nucleic Acids Res. 2005 Jan 1;33 Database Issue:D562-6"
    reference_series
    [1] "GSE781"
    sample_count
    [1] "17"
    sample_id
     [1] "GSM11815,GSM11832,GSM12069,GSM12083,GSM12101,GSM12106,GSM12274,GSM12299,GSM12412"
     [2] "GSM11810,GSM11827,GSM12078,GSM12099,GSM12269,GSM12287,GSM12301,GSM12448"
     [3] "GSM11810,GSM11815"
     [4] "GSM11827,GSM11832"
     [5] "GSM12069,GSM12078"
     [6] "GSM12083,GSM12099"
     [7] "GSM12101"
     [8] "GSM12106"
     [9] "GSM12269"
    [10] "GSM12274,GSM12287"
    [11] "GSM12299,GSM12301"
    [12] "GSM12412,GSM12448"
    sample_organism
    [1] "Homo sapiens"
    sample_type
    [1] "RNA"
    title
    [1] "Renal clear cell carcinoma (HG-U133B)"
    type
     [1] "Expression profiling by array" "disease state"
     [3] "disease state"                 "individual"
     [5] "individual"                    "individual"
     [7] "individual"                    "individual"
     [9] "individual"                    "individual"
    [11] "individual"                    "individual"
    [13] "individual"
    update_date
    [1] "Mar 04 2004"
    value_type
    [1] "count"
    web_link
    [1] "http://www.ncbi.nlm.nih.gov/geo"
    An object of class "GEODataTable"
    ****** Column Descriptions ******
         sample disease.state individual
    1  GSM11815           RCC        035
    2  GSM11832           RCC        023
    3  GSM12069           RCC        001
    4  GSM12083           RCC        005
    5  GSM12101           RCC        011
    6  GSM12106           RCC        032
    7  GSM12274           RCC          2
    8  GSM12299           RCC          3
    9  GSM12412           RCC          4
    10 GSM11810        normal        035
    11 GSM11827        normal        023
    12 GSM12078        normal        001
    13 GSM12099        normal        005
    14 GSM12269        normal          1
    15 GSM12287        normal          2
    16 GSM12301        normal          3
    17 GSM12448        normal          4
                                                                                                                                           description
    1             Value for GSM11815: C035 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    2             Value for GSM11832: C023 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    3             Value for GSM12069: C001 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    4             Value for GSM12083: C005 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    5             Value for GSM12101: C011 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    6             Value for GSM12106: C032 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    7               Value for GSM12274: C2 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    8               Value for GSM12299: C3 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    9               Value for GSM12412: C4 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
    10      Value for GSM11810: N035 Normal Human Kidney U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    11      Value for GSM11827: N023 Normal Human Kidney U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    12      Value for GSM12078: N001 Normal Human Kidney U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    13      Value for GSM12099: N005 Normal Human Kidney U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    14        Value for GSM12269: N1 Normal Human Kidney U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    15 Value for GSM12287: N2 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    16 Value for GSM12301: N3 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    17 Value for GSM12448: N4 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from normal tissue adjacent to Renal Cell Carcinoma
    ****** Data Table ******

GDS objects can be converted to Bioconductor data structures:

``` r
# Convert to ExpressionSet (with log2 transformation)
eset <- GDS2eSet(gds, do.log2=TRUE)
eset
```

    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 22645 features, 17 samples
      element names: exprs
    protocolData: none
    phenoData
      sampleNames: GSM11815 GSM11832 ... GSM12448 (17 total)
      varLabels: sample disease.state individual description
      varMetadata: labelDescription
    featureData
      featureNames: 200000_s_at 200001_at ... AFFX-TrpnX-M_at (22645 total)
      fvarLabels: ID Gene title ... GO:Component ID (21 total)
      fvarMetadata: Column labelDescription
    experimentData: use 'experimentData(object)'
      pubMedIds: 14641932
    Annotation:  

``` r
# Or to a limma MAList
malist <- GDS2MA(gds)
class(malist)
```

    [1] "MAList"
    attr(,"package")
    [1] "limma"

These conversions are particularly useful for integrating older GEO
datasets into modern analytical workflows.

## Advanced Features

### Getting GSE Data Tables

Some GSE records contain data tables with important metadata not
captured in the standard GSE structure:

``` r
# Get data tables from GSE3494
dt_list <- getGSEDataTables("GSE3494")
```

    Rows: 251 Columns: 12
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (6): X1, X2, X5, X6, X7, X10
    dbl (6): X3, X4, X8, X9, X11, X12

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    Rows: 502 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (3): X1, X2, X3

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
names(dt_list)
```

    NULL

``` r
head(dt_list[[1]])
```

    # A tibble: 6 × 12
      `INDEX (ID)` p53 seq mut status (p53+=mutant; p53-=wt…¹ p53 DLDA classifier …²
      <chr>        <chr>                                                       <dbl>
    1 X101B88      p53+                                                            1
    2 X102B06      p53+                                                            1
    3 X104B91      p53+                                                            0
    4 X110B34      p53+                                                            1
    5 X111B51      p53+                                                            1
    6 X127B00      p53+                                                            1
    # ℹ abbreviated names: ¹​`p53 seq mut status (p53+=mutant; p53-=wt)`,
    #   ²​`p53 DLDA classifier result (0=wt-like, 1=mt-like)`
    # ℹ 9 more variables: `DLDA error (1=yes, 0=no)` <dbl>,
    #   `Elston histologic grade` <chr>, `ER status` <chr>, `PgR status` <chr>,
    #   `age at diagnosis` <dbl>, `tumor size (mm)` <dbl>,
    #   `Lymph node status` <chr>,
    #   `DSS TIME (Disease-Specific Survival Time in years)` <dbl>, …

### Working with GPL Platforms

Platform records (GPL) contain important probe annotations:

``` r
# Get a platform record
gpl <- getGEO("GPL96")
head(Table(gpl)[, 1:5])
```

             ID GB_ACC SPOT_ID Species Scientific Name Annotation Date
    1 1007_s_at U48705                    Homo sapiens     Oct 6, 2014
    2   1053_at M87338                    Homo sapiens     Oct 6, 2014
    3    117_at X51757                    Homo sapiens     Oct 6, 2014
    4    121_at X69699                    Homo sapiens     Oct 6, 2014
    5 1255_g_at L36861                    Homo sapiens     Oct 6, 2014
    6   1294_at L13852                    Homo sapiens     Oct 6, 2014

When retrieving GSE records, GEOquery can automatically include GPL
annotation:

``` r
# Get GSE with GPL annotation
gse_with_gpl <- getGEO("GSE2553", AnnotGPL=TRUE)
head(fData(gse_with_gpl[[1]]))
```

## Reporting Bugs and Contributing

As GEO continues to evolve, GEOquery adapts to support new features and
data types. If you encounter issues:

1.  Check the [Bioconductor Support
    site](https://support.bioconductor.org/)
2.  Report bugs on [GitHub](https://github.com/seandavi/GEOquery/issues)
3.  Consider contributing via pull requests

## Citing GEOquery

If you use GEOquery in your research, please cite:

``` r
citation("GEOquery")
```

    Please cite the following if utilizing the GEOquery software:

      Davis S, Meltzer P (2007). "GEOquery: a bridge between the Gene
      Expression Omnibus (GEO) and BioConductor." _Bioinformatics_, *14*,
      1846-1847. doi:10.1093/bioinformatics/btm254
      <https://doi.org/10.1093/bioinformatics/btm254>.

    A BibTeX entry for LaTeX users is

      @Article{,
        author = {Sean Davis and Paul Meltzer},
        title = {GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor},
        journal = {Bioinformatics},
        year = {2007},
        volume = {14},
        pages = {1846--1847},
        doi = {10.1093/bioinformatics/btm254},
      }

## Session Information

``` r
sessionInfo()
```

    R version 4.5.2 (2025-10-31)
    Platform: aarch64-apple-darwin20
    Running under: macOS Sequoia 15.7.1

    Matrix products: default
    BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: UTC
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] knitr_1.50          GEOquery_2.77.6     Biobase_2.70.0
    [4] BiocGenerics_0.56.0 generics_0.1.4

    loaded via a namespace (and not attached):
     [1] SummarizedExperiment_1.40.0 xfun_0.54
     [3] httr2_1.2.1                 lattice_0.22-7
     [5] tzdb_0.5.0                  vctrs_0.6.5
     [7] tools_4.5.2                 parallel_4.5.2
     [9] stats4_4.5.2                curl_7.0.0
    [11] tibble_3.3.0                pkgconfig_2.0.3
    [13] R.oo_1.27.1                 Matrix_1.7-4
    [15] data.table_1.17.8           rentrez_1.2.4
    [17] S4Vectors_0.48.0            lifecycle_1.0.4
    [19] compiler_4.5.2              stringr_1.6.0
    [21] statmod_1.5.1               Seqinfo_1.0.0
    [23] htmltools_0.5.8.1           yaml_2.3.10
    [25] pillar_1.11.1               crayon_1.5.3
    [27] tidyr_1.3.1                 R.utils_2.13.0
    [29] DelayedArray_0.36.0         limma_3.66.0
    [31] abind_1.4-8                 tidyselect_1.2.1
    [33] rvest_1.0.5                 digest_0.6.39
    [35] stringi_1.8.7               dplyr_1.1.4
    [37] purrr_1.2.0                 fastmap_1.2.0
    [39] grid_4.5.2                  cli_3.6.5
    [41] SparseArray_1.10.2          magrittr_2.0.4
    [43] S4Arrays_1.10.0             utf8_1.2.6
    [45] XML_3.99-0.20               readr_2.1.6
    [47] withr_3.0.2                 rappdirs_0.3.3
    [49] bit64_4.6.0-1               rmarkdown_2.30
    [51] XVector_0.50.0              httr_1.4.7
    [53] matrixStats_1.5.0           bit_4.6.0
    [55] R.methodsS3_1.8.2           hms_1.1.4
    [57] evaluate_1.0.5              GenomicRanges_1.62.0
    [59] IRanges_2.44.0              rlang_1.1.6
    [61] glue_1.8.0                  selectr_0.5-0
    [63] xml2_1.5.0                  vroom_1.6.6
    [65] jsonlite_2.0.0              R6_2.6.1
    [67] MatrixGenerics_1.22.0      
