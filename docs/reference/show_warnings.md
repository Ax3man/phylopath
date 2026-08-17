# Print out warnings from a phylopath analysis.

Use this function after running
[`phylo_path()`](https://ax3man.github.io/phylopath/reference/phylo_path.md)
to conveniently print any generated warnings to the screen. You can
either provide no arguments, which will only work if you run it directly
after the analysis, or you have to provide the `phylopath` object
manually.

## Usage

``` r
show_warnings(phylopath = NULL)
```

## Arguments

- phylopath:

  A `phylopath` object of which the warnings should be printed.
