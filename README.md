# phylopath

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/phylopath)](https://cran.r-project.org/package=phylopath) [![Cran downloads](http://cranlogs.r-pkg.org/badges/grand-total/phylopath)](http://cran.rstudio.com/web/packages/phylopath/index.html) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

This package implements phylogenetic path analysis in R.

Install the package using:

```{r}
install.packages("phylopath")
```
You may need to install the downstream dependency `graph` from Bioconductor:

```{r}
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("graph")
```

Or if you'd like to install the development version (here on github), use:

```{r}
remotes::install_github("Ax3man/phylopath")
```

It's easiest to start on the [website](https://ax3man.github.io/phylopath) and first read the introduction [here](https://ax3man.github.io/phylopath/articles/intro_to_phylopath.html), or read the [paper in PeerJ](https://doi.org/10.7717/peerj.4718).

If you find any problems, or if you have suggestions for improvements, please file those under [issues](/issue). PRs welcome.
