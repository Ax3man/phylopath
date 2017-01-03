phylopath 0.2.2
--------------------------------------------------------------------------------

* Fewer models are now fitted when using `phylo_path()`, since any duplicated
  independence statements are now only fitted once. This leads to a significant
  reduction in running time in many cases, especially when many models are
  conscidered.
  
* Implemented support for parallel processing in `phylo_path()` using the
  `parallel` argument.
  
* `phylo_path()` now shows a progress bar. 

* New function added (`choice()`) that is a very simple wrapper around 
  `est_DAG()`. It adds to `best()` and `average()` by allowing for choosing
  any model as the final model, and encourages users to not always pick the 
  lowest CICc model.
  
* Prepared plotting functions for new release of `DiagrammeR`, v0.9 now
  required.


phylopath 0.2.1
--------------------------------------------------------------------------------

* IMPORTANT: Faulty model averaging has been fixed. This was often introduced
  due to differences in matrix ordering. Avering results from versions before
  0.2.1 should NOT be trusted.

* Using `ape::corBrownian()` no longer returns an error.

* Averaging is less likely to fail due to errors in `nlme::intervals()`.


phylopath 0.2.0
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
  