## Test environments
* local Windows 10 install, R 3.3.3
* Development version using winbuilder r72375

## R CMD check results
0 errors | 0 warnings | 0 notes

## CRAN check:

There is currently 1 NOTE on some builds:

Result: NOTE 
    Namespace in Imports field not imported from: ‘ape’
     All declared Imports should be used. 

While no function from ape is directly used, a function from ape is used as a
default argument in the phylo_path function.

## Downstream dependencies
The package does not have downstream dependencies.
