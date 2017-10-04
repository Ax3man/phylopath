phylopath 0.3.1.9000
--------------------------------------------------------------------------------

* `phylo_path` and all related methods now deal automatically with both 
  continuous and binary data. All separate binary functions and methods have
  disappeared as they are no longer needed.

* Added `plot()` method for `phylopath.summary` objects, that shows the weights
  and p-values for the different models.

* Added a `define_model_set()` convenience function for building models, that 
  avoids repeated calls to `DAG()` and has an argument to supply paths that are 
  shared between all your models. It is not needed to specify isolate variables.
  Old code using `DAG()` continues to work as normal.

* Added support for additional arguments passed to `gls` from `phylo_path`. This
  can be helpful, for example, for setting the fitting method to maximum 
  likelihood (`method = "ML"`).

phylopath 0.3.1
--------------------------------------------------------------------------------

####Bugfixes:

* The package broke due to an update of `purrr`, but has now been fixed 
  (reported by Christoph Liedtke, @hcliedtke).
  
* The package depends on a recent version of `nlme`, but this wasn't specified.
  All package versions of dependencies are now defined. (reported by @ManuelaGonzalez)

phylopath 0.3.0
--------------------------------------------------------------------------------

* Added support for completely binary models, that are fitted with 
  `ape::binaryPGLMM`. Use `phylo_path_binary()` to compare models. `average()`,
  `best()` and `choice()` are now S3 generics and will handle both continuous
  and binary versions. Usage is designed to be as close to the continuous version
  as possible. `est_DAG_binary()` powers the binary S3 methods.

* All plot functions that used `DiagrammeR` now use `ggraph` instead. This gives
  much more control over the positioning of the nodes, and allows to plot 
  multiple models at once. Exporting plots also becomes much easier.

* You can now plot a list of causal models with `plot_model_set()`. This 
  creates a faceted plot where all nodes are kept in the same location, which 
  makes it easier to spot how models are different.

phylopath 0.2.3
--------------------------------------------------------------------------------

* If there are any `NA` values in `data` for the variables in `models`, these
  rows are now dropped from `data` with a message. Use `na.rm = FALSE` to revert
  to the old behavior.

* When PGLS models fail, an informative error is now returned to the user.

* `phylo_path()` now checks for row.names that line up with the tree tip labels.
  If the tree contains surplus species, it gets pruned to size with a message.
  This includes cases where species are dropped due to missing values.
  
* `citation()` now correctly refers to the methods paper by Von Hardenberg &
  Gonzalez-Voyer first and the package second.

phylopath 0.2.2
--------------------------------------------------------------------------------

* Fewer models are now fitted when using `phylo_path()`, since any duplicated
  independence statements are now only fitted once. This leads to a significant
  reduction in running time in many cases, especially when many models are
  considered.
  
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
  due to differences in matrix ordering. Averaging results from versions before
  0.2.1 should NOT be trusted.

* Using `ape::corBrownian()` no longer returns an error.

* Averaging is less likely to fail due to errors in `nlme::intervals()`.


phylopath 0.2.0
--------------------------------------------------------------------------------

* `phylo_path()` has become more streamlined with functionality moved to other
  functions. The `phylopath` object now contains all necessary models and data,
  `summary()` is used to obtain the results table, and `best()` and `average()` 
  are used to extract and fit the best or average model. See the vignette for
  details.

* Model averaging for arbitrary models is now possible with `average_DAGs()`.

* Model averaging now supports both conditional and full model averaging.

* Both the old `est_DAG()` and the new `average_DAGs()` now return objects of a
  new class `fitted_DAG`, that has it's separate `plot` method. The `plot` 
  method for objects of class `DAG` has been simplified.

* Model averaging now returns standard errors and confidence intervals based on
  the `MuMIn` package (issue #1).
  
* A new function `plot_coefs` for plotting regression coefficients and their
  confidence intervals has been added.
  
