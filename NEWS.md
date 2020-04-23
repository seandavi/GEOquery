CHANGES IN VERSION 2.55.1
---------------------------

Bug fixes:

    * Fix test failure due to stricter boolean check lengths [#568ae94fb95, from @russHyde]

    * Fix failure on empty file [#88741521a0e, from @vlakam]

CHANGES IN VERSION 2.51.1
---------------------------

Improvements:

    * Add formal experimentData slot to GSE records. Of class MIAME. [from @vlakam]
      
CHANGES IN VERSION 2.47.1
------------------------------

Bug fixes:

    * Fixes problems with intermittent connection issues (which were not intermittend connection problems, it seems)

CHANGES IN VERSION 2.45.2
-----------------------------

Improvements:

    * GPL parsing 4-5x faster
    * GSM parsing 3x faster
    * GSEMatrix parsing is much smarter with respect to sample characteristics. In short, for GSEs where sample characteristics are actually used, the pData should have nice, neat column headers with the phenodata keys and values in the columns, including correct handling of missing values, etc.

CHANGES IN VERSION 2.45.1
-----------------------------

Bug fixes

    * getDirectoryListing fixed to deal with changes
      to NCBI server listing formats

CHANGES IN VERSION 2.45
---------------------------

Improvements:

    * GDS parsing is 2-3x faster
    * GSEMatrix parsing is 2-3x faster

CHANGES IN VERSION 2.36
------------------------
    
New Features

    * New, faster SOFT format parsing (Leonardo Gama)
    * Turned on unit tests in Travis CI
    * Test coverage metrics added
    
Bug fixes

    * default download method no longer assumes that curl is installed on linux
    * GSEMatrix parsing from file now finds cached GPLs
