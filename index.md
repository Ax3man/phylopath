# Easy phylogenetic path analysis in R

Use the `phylopath` package for an easy to use framework to perform
phylogenetic path analysis (PPA).

PPA can be used to compare support for competing causal models of trait
evolution, while taking shared ancestry into account. All you need
is: 1. A clear set of models to test. 2. A data set of species with
trait values. 3. A phylogeny of your species.

For a complete worked example, click “Get Started” above, or see the
[PeerJ paper](https://doi.org/10.7717/peerj.4718).

## Installation

    install.packages("BiocManager")
    BiocManager::install("phylopath")

`phylopath` relies on `ggm`, which in turn needs `graph` from
Bioconductor. `BiocManager::install()` installs from CRAN as well as
Bioconductor, so this single command obtains all three. A plain
`install.packages("phylopath")` works only if you already have `graph`.

This method was developed by Von Hardenberg and Gonzalez-Voyer. See
`citation()` for info on correct citations.

<img src="man/figures/unnamed-chunk-2-1.png" alt="A grid of nine small causal diagrams, one per candidate model, each showing directed arrows between the five rhinograde traits body mass, litter size, nose length, dry days and range size." width="600px" style="display: block; margin: auto;" />

<img src="man/figures/unnamed-chunk-3-1.png" alt="The fitted best supported causal model, with arrows labelled by their standardized path coefficients and drawn with a width proportional to the strength of the effect." width="600px" style="display: block; margin: auto;" />

`phylopath` has been used by [&gt;200
publications](https://scholar.google.ca/scholar?oi=bibs&hl=en&cites=5933615079034924484)
so far!
