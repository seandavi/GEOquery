# GSE Supplemental file listing

The GEO Series records often have one or more supplemental files. In
most cases, those files are archived as '.tar' files, the contents of
which are only available in a file listing file not present on the
website for download.

## Usage

``` r
getGEOSeriesFileListing(GSE)
```

## Arguments

- GSE:

  character(1) the GSE accession

## Value

A data.frame with 5 columns. See example.

## Details

This function reads that file listing file and returns the results as a
data.frame.

## Examples

``` r
getGEOSeriesFileListing('GSE288770')
#> Rows: 17 Columns: 5
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (4): #Archive/File, Name, Time, Type
#> dbl (1): Size
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> # A tibble: 17 × 5
#>    archive_or_file name                                       time    size type 
#>    <chr>           <chr>                                      <chr>  <dbl> <chr>
#>  1 Archive         GSE288770_RAW.tar                          08/0… 1.20e8 TAR  
#>  2 File            GSM8775062_vitro_contol_filtered_feature_… 02/0… 3.01e7 H5   
#>  3 File            GSM8775063_vitro_vasc_filtered_feature_bc… 02/0… 1.54e7 H5   
#>  4 File            GSM8775064_vitro_neuro_filtered_feature_b… 02/0… 2.72e7 H5   
#>  5 File            GSM8775065_vitro_odont_filtered_feature_b… 02/0… 2.22e7 H5   
#>  6 File            GSM8775066_vivo_day7_1_barcodes.tsv.gz     02/0… 3.36e3 TSV  
#>  7 File            GSM8775066_vivo_day7_1_features.tsv.gz     02/0… 3.34e5 TSV  
#>  8 File            GSM8775066_vivo_day7_1_matrix.mtx.gz       02/0… 2.98e6 MTX  
#>  9 File            GSM8775067_vivo_day7_2_barcodes.tsv.gz     02/0… 7.95e3 TSV  
#> 10 File            GSM8775067_vivo_day7_2_features.tsv.gz     02/0… 3.34e5 TSV  
#> 11 File            GSM8775067_vivo_day7_2_matrix.mtx.gz       02/0… 1.28e7 MTX  
#> 12 File            GSM8775068_vivo_day21_1_barcodes.tsv.gz    02/0… 9.71e3 TSV  
#> 13 File            GSM8775068_vivo_day21_1_features.tsv.gz    02/0… 3.34e5 TSV  
#> 14 File            GSM8775068_vivo_day21_1_matrix.mtx.gz      02/0… 5.72e6 MTX  
#> 15 File            GSM8775069_vivo_day21_2_barcodes.tsv.gz    02/0… 2.1 e3 TSV  
#> 16 File            GSM8775069_vivo_day21_2_features.tsv.gz    02/0… 3.34e5 TSV  
#> 17 File            GSM8775069_vivo_day21_2_matrix.mtx.gz      02/0… 2.41e6 MTX  
```
