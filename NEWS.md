pylopath 0.0.2
--------------------------------------------------------------------------------

* `phylo_path()` has become more streamlined with functionality moved to other
  functions. The `phylopath` object now containts all necessary models and data,
  `summary()` is used to obtain the results table, and `best()` and `average()` 
  are used to extract and fit the best or average model. See the vignette for
  details.

* Model averaging for arbitrary models is now possible with `average_DAGs()`.

* Model averaging now supports both conditional and full model averaging.

* Both the old `est_DAG()` and the new `average_DAGs()` now return objects of a
  new class `fitted_DAG`, that has it's seperate `plot` method. The `plot` 
  method for objects of class `DAG` has been simplified.

* Model averaging now returns standard errors and confidence intervals based on
  the `MuMIn` package (issue #1).
  
* A new function `plot_coefs` for plotting regression coefficients and their
  confidence intervals has been added.