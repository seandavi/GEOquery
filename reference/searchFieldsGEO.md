# Provide a list of possible search fields for GEO search

Provide a list of possible search fields for GEO search

## Usage

``` r
searchFieldsGEO()
```

## Value

a data.frame with names of possible search fields for GEO search as well
as descriptions, data types, etc. for each field. Fields are in rows and
their properties are in columns.

## See also

[`searchGEO`](http://seandavi.github.io/GEOquery/reference/searchGEO.md)

## Examples

``` r
searchFieldsGEO()
#>    Name     FullName  Description TermCount IsDate IsNumerical SingleToken
#> 1   ALL   All Fields All term....  45152635      N           N           N
#> 2   UID          UID Unique n....         0      N           Y           Y
#> 3  FILT       Filter Limits t....        71      N           N           Y
#> 4  ORGN     Organism exploded....     75623      N           N           Y
#> 5  ACCN GEO Acce.... accessio....  19303976      N           N           Y
#> 6  TITL        Title Words in....   9946487      N           N           Y
#> 7  DESC  Description Text fro....  10548780      N           N           Y
#> 8  SFIL Suppleme.... Suppleme....       255      N           N           Y
#> 9  ETYP   Entry Type Entry ty....         4      N           N           Y
#> 10 STYP  Sample Type  Sample type         9      N           N           Y
#> 11 VTYP Sample V.... type of ....         7      N           N           Y
#> 12 PTYP Platform.... Platform....        17      N           N           Y
#> 13 GTYP DataSet Type type of ....        27      N           N           Y
#> 14 NSAM Number o.... Number o....      2134      N           Y           Y
#> 15  SRC Sample S.... sample s....    461250      N           N           Y
#> 16 AUTH       Author author o....   1300542      N           N           Y
#> 17 INST Submitte.... institut....     25314      N           N           Y
#> 18 NPRO Number o.... number o....      7258      N           Y           Y
#> 19 SSTP Subset V.... subset v....        24      N           N           Y
#> 20 SSDE Subset D.... subset d....      7535      N           N           Y
#> 21 GEID Reporter.... name or ....   2840498      N           N           Y
#> 22 PDAT Publicat.... publicat....      8453      Y           N           Y
#> 23 UDAT  Update Date         date      7570      Y           N           Y
#> 24 TAGL   Tag Length Tag/Sign....         9      N           N           Y
#> 25 RGSE Related .... Related ....     30135      N           N           Y
#> 26 RGPL Related .... Related ....    271071      N           N           Y
#> 27 MESH   MeSH Terms Medical ....     17643      N           N           Y
#> 28 PROJ      Project      Project        10      N           N           Y
#> 29 ATNM Attribut.... Attribut....     49222      N           N           Y
#> 30 ATTR    Attribute    Attribute   2763893      N           N           Y
#> 31 PROP   Properties   Properties         3      N           N           Y
#>    Hierarchy IsHidden
#> 1          N        N
#> 2          N        Y
#> 3          N        N
#> 4          Y        N
#> 5          N        N
#> 6          N        N
#> 7          N        N
#> 8          N        N
#> 9          N        N
#> 10         N        N
#> 11         N        N
#> 12         N        N
#> 13         N        N
#> 14         N        N
#> 15         N        N
#> 16         N        N
#> 17         N        N
#> 18         N        N
#> 19         N        N
#> 20         N        N
#> 21         N        N
#> 22         N        N
#> 23         N        N
#> 24         N        N
#> 25         N        N
#> 26         N        N
#> 27         Y        N
#> 28         N        N
#> 29         N        N
#> 30         N        N
#> 31         N        N
```
