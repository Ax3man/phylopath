## Test environments
* local macOS 26.5.1 install, R 4.6.1
* GitHub Actions: macOS (release), Windows (release, devel), Ubuntu (devel, release, oldrel-1)

## Local R CMD check results
0 errors | 0 warnings | 0 notes

## Wincheck results
1 NOTE generated:
    Possibly misspelled words in DESCRIPTION:
        BiocManager (13:60, 14:5)
        phylopath (14:27)
These words are not misspelled.

## CRAN Package Check Results
This release fixes the current NOTE on the r-devel flavours about deprecated
arguments to `structure()` in the test code ("'.Dim' should be changed to 'dim'",
"'.Dimnames' should be changed to 'dimnames'"). Those calls no longer occur
anywhere in the packages.

## Downstream dependencies
The package has 2 downstream dependencies, which were checked with revdepcheck and passed.
