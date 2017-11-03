## Test environments
* local Windows 10 install (x86_64-w64-mingw32/x64), R 3.4.2
* Development version using winbuilder r73242

## Local R CMD check results
0 errors | 0 warnings | 0 notes

## winbuilder Dev results
0 errors | 0 warnings | 0 notes

## CRAN check:

There are currently errors on r-oldrel builds because of a version dependency for the parallel
package. I've removed the specific version dependency, so these errors should dissapear with this
release.

## Downstream dependencies
The package does not have downstream dependencies.
