## Test environments
* local Windows 10 install (x86_64-w64-mingw32/x64), R 3.4.1
* Development version using winbuilder r73068

## Local R CMD check results
0 errors | 0 warnings | 0 notes

## winbuilder Dev results
0 errors | 0 warnings | 1 note

The note is about possible misspellings. The spelling is correct.

## CRAN check:

There is currently errors on most builds. This was caused by an update of `purrr`, of which I was
not informed. This happened because both `purrr` and `phylopath` pushed updates around the same
time and the maintainers of `purrr` did not see the newly found errors.

The current release fixes these errors.

## Downstream dependencies
The package does not have downstream dependencies.
