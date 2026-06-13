# Get GSE data tables from GEO into R data structures.

In some cases, instead of individual sample records (GSM) containing
information regarding sample phenotypes, the GEO Series contains that
information in an attached data table. And example is given by GSE3494
where there are two data tables with important information contained
within them. Using getGEO with the standard parameters downloads the
GSEMatrix file which, unfortunately, does not contain the information in
the data tables. This function simply downloads the “header” information
from the GSE record and parses out the data tables into R data.frames.

## Usage

``` r
getGSEDataTables(GSE)
```

## Arguments

- GSE:

  The GSE identifier, such as “GSE3494”.

## Value

A list of data.frames.

## See also

[`getGEO`](http://seandavi.github.io/GEOquery/reference/getGEO.md)

## Author

Sean Davis <sdavis2@mail.nih.gov>

## Examples

``` r

dfl = getGSEDataTables('GSE3494')
#> Rows: 251 Columns: 12
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (6): X1, X2, X5, X6, X7, X10
#> dbl (6): X3, X4, X8, X9, X11, X12
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 502 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (3): X1, X2, X3
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
lapply(dfl,head)
#> [[1]]
#> # A tibble: 6 × 12
#>   `INDEX (ID)` p53 seq mut status (p53+=mutant; p53-=wt…¹ p53 DLDA classifier …²
#>   <chr>        <chr>                                                       <dbl>
#> 1 X101B88      p53+                                                            1
#> 2 X102B06      p53+                                                            1
#> 3 X104B91      p53+                                                            0
#> 4 X110B34      p53+                                                            1
#> 5 X111B51      p53+                                                            1
#> 6 X127B00      p53+                                                            1
#> # ℹ abbreviated names: ¹​`p53 seq mut status (p53+=mutant; p53-=wt)`,
#> #   ²​`p53 DLDA classifier result (0=wt-like, 1=mt-like)`
#> # ℹ 9 more variables: `DLDA error (1=yes, 0=no)` <dbl>,
#> #   `Elston histologic grade` <chr>, `ER status` <chr>, `PgR status` <chr>,
#> #   `age at diagnosis` <dbl>, `tumor size (mm)` <dbl>,
#> #   `Lymph node status` <chr>,
#> #   `DSS TIME (Disease-Specific Survival Time in years)` <dbl>, …
#> 
#> [[2]]
#> # A tibble: 6 × 3
#>   `GEO Sample Accession #` `Patient ID` `Affy platform`
#>   <chr>                    <chr>        <chr>          
#> 1 GSM79114                 X100B08      HG-U133A       
#> 2 GSM79115                 X101B88      HG-U133A       
#> 3 GSM79116                 X102B06      HG-U133A       
#> 4 GSM79117                 X103B41      HG-U133A       
#> 5 GSM79118                 X104B91      HG-U133A       
#> 6 GSM79119                 X105B13      HG-U133A       
#> 

```
