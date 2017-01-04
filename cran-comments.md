## Test environments
* local Windows 10 install, R 3.3.2

## R CMD check results
0 errors | 0 warnings | 0 notes

## CRAN check:

There are currently 1 ERROR and 1 NOTE:

Result: ERROR 
    
  Running examples in ‘phylopath-Ex.R’ failed

This ERROR occurs due to an update in the DiagrammeR package to v0.9.0. This
release fixes all errors.

Result: NOTE 
    Namespace in Imports field not imported from: ‘ape’
     All declared Imports should be used. 

While no function from ape is directly used, a function from ape is used as a
default argument in the phylo_path function.

## Downstream dependencies
The package does not have downstream dependencies.
