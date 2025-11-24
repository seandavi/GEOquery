# Convert a GDS data structure to a BioConductor data structure

Functions to take a GDS data structure from getGEO and coerce it to
limma MALists or ExpressionSets.

## Arguments

- GDS:

  The GDS datastructure returned by getGEO

- do.log2:

  Boolean, should the data in the GDS be log2 transformed before
  inserting into the new data structure

- GPL:

  Either a GPL data structure (from a call to getGEO) or NULL. If NULL,
  this will cause a call to getGEO to produce a GPL. The gene
  information from the GPL is then used to construct the `genes` slot of
  the resulting limma `MAList` object or the `featureData` slot of the
  `ExpressionSet` instance.

- AnnotGPL:

  In general, the annotation GPL files will be available for GDS
  records, so the default is to use these files over the user-submitted
  GPL files

- getGPL:

  A boolean defaulting to TRUE as to whether or not to download and
  include GPL information when converting to ExpressionSet or MAList.
  You may want to set this to FALSE if you know that you are going to
  annotate your featureData using Bioconductor tools rather than relying
  on information provided through NCBI GEO. Download times can also be
  greatly reduced by specifying FALSE.

## Value

- GDS2MA:

  A limma MAList

- GDS2eSet:

  An ExpressionSet object

## Details

This function just rearranges one data structure into another. For GDS,
it also deals appropriately with making the 'targets' list item for the
limma data structure and the phenoData slot of ExpressionSets.

## References

See the limma and ExpressionSet help in the appropriate packages

## Author

Sean Davis

## Examples

``` r

if (FALSE) gds505 <- getGEO('GDS505') # \dontrun{}
if (FALSE) MA <- GDS2MA(gds505) # \dontrun{}
if (FALSE) eset <- GDS2eSet(gds505) # \dontrun{}

```
