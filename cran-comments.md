## Test environments
* local Windows 10 install (x86_64-w64-mingw32/x64), R 3.4.1
* Development version using winbuilder r72891

## Local R CMD check results
0 errors | 0 warnings | 0 notes

## winbuilder Dev results
0 errors | 0 warnings | 1 note

The note about possible misspellings. The spelling is correct.

## CRAN check:

There is currently 1 NOTE on some builds:

Result: NOTE 
    Namespace in Imports field not imported from: ‘ape’
     All declared Imports should be used. 

ape is now directly used, so this release fixes that NOTE.

## Downstream dependencies
The package does not have downstream dependencies.
