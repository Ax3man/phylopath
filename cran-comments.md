## Test environments
* local Windows 10 install, R 3.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * checking CRAN incoming feasibility ... NOTE
Maintainer: 'Wouter van der Bijl <wouter.van.der.bijl@zoologi.su.se>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Hardenberg (8:52)
  Phylogenetic (3:16)
  Von (8:48)
  Voyer (8:76)
  phylogenetic (8:5)

These words are not mis-spelled, but include names and 'phylogenetic' analysis
refers to analyses that control for interdependence of species due to the tree
of life.

## Downstream dependencies
The package does not have downstream dependencies.

## Comments regarding previous rejected submission
Thanks! I've have made these changes:
* doi has been added to both DESCRIPTION as well as data documentation.
* I have added examples to all exported functions with documentation.
