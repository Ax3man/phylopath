# Directed acyclic graphs (DAGs)

This function is a simple wrapper around the function from the `ggm`
package with the same name. The only differences are that the `order`
argument defaults to `TRUE` and that it adds a `DAG` class for easy
plotting. Typically, one would use
[`define_model_set()`](https://ax3man.github.io/phylopath/reference/define_model_set.md)
to create models for use with the `phylopath` package.

## Usage

``` r
DAG(..., order = TRUE)
```

## Arguments

- ...:

  a sequence of model formulae

- order:

  logical, defaulting to `TRUE`. If `TRUE` the nodes of the DAG are
  permuted according to the topological order. If `FALSE` the nodes are
  in the order they first appear in the model formulae (from left to
  right). For use in the `phylopath` package, this should always be kept
  to `TRUE`, but the argument is available to avoid potential problems
  with masking the function from other packages.

## Value

An object of classes `matrix` and `DAG`

## Details

Supply a formulas for the model as arguments. Formulas should be of the
form
``` child ~ parent`` and describe each path in your model. Multiple children of a single parent can be combined into a single formula:  ```child
~ parent1 +
parent2`. Finally, an isolate (unconnected variable) can be included as being connected to itself: `isolate
~ isolate\`.

## Examples

``` r
  # Use formula notation to create DAGs:
  plot(DAG(A~B, B~C))

  # Use + to easily add multiple parents to a node:
  plot(DAG(A~B+C))

  # Add a node as it's own parent to create an isolate:
  plot(DAG(A~B+C, D~D))
```
