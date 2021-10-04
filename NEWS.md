phylopath 1.1.3
--------------------------------------------------------------------------------

* Fixed a bug that made `phylo_path` fail to pass additional (...) 
  arguments correctly to `phylolm`.

* Add informative error when trying to plot a DAG without any paths.

* Updated plotting functions to work with new `ggraph` releases.

* Fixed regression with parallel usage of `phylo_path` due to an S3 
  inheritance issue on the cluster (#16, thanks Simon Greenhill for the report).

phylopath 1.1.2
--------------------------------------------------------------------------------

* Prepare for R v4.0.0.

* Bug fix: Very low p-values could cause underflow and result in infinite C 
  statistics. All p-values are now set to be at least the size of the machine
  accuracy (i.e. 2 * 10^-16).
  
* Warnings are now again correctly reported.

phylopath 1.1.1
--------------------------------------------------------------------------------

* Prepared for next release of `ggraph`.

phylopath 1.1.0
--------------------------------------------------------------------------------

* Bug fix: It was possible to get CICc values in the summary output that were
  not valid. Specifically, to calculate CICc there is a division by 
  `(n - 1 - q)`, where `n` is the number of observations (species) and `q` the 
  number of parameters in the causal model. This could lead to infinite CICc
  when `n == q + 1`, or a flipped of CICc when `n < q + 1`. This would typically
  only occur when attempting to fit models with very few species (e.g. < 10).
  
  New behavior is to set CICc to `NA` when `n` is insufficient, and to give a
  warning.
  
* Removed dependencies `dplyr` and `tidyr`, but added `tibble`.

phylopath 1.0.2
--------------------------------------------------------------------------------

* Prepare for dplyr 0.8.0 release.

phylopath 1.0.1
--------------------------------------------------------------------------------

* Fixed bug that would return the wrong model in some error messages.

* Improved reporting of warnings, and a `show_warnings()` function has been 
  added.

* Citation info now points to the [PeerJ paper](https://doi.org/10.7717/peerj.4718).

phylopath 1.0.0
--------------------------------------------------------------------------------

* Citation info now points to the [bioRxiv paper](https://www.biorxiv.org/content/10.1101/212068v1).

* All modeling functions now completely rely on the `phylolm` package, and no 
  longer use `ape`. This is a major change, that will possibly change the 
  outcomes of some of your existing analyses (as can happen when chaning 
  the modeling package). There are, however, several good reasons to make this 
  change, which I think make it worth the trouble. Firstly, the package is much
  faster for large trees, and this effect is compounded in `phylopath` because 
  one may have to fit a few dozen models. Secondly, I think it is important to 
  have confidence intervals around the regression coefficients, and those were 
  not available for `ape::binaryPGLMM`. Thirdly, `phylolm` makes it easy to use 
  a larger variety of models of evolution, including two versions of OU and 
  early burst, which can be simply set using the `model` parameter. Lastly, the
  `phylolm()` and `phyloglm()` functions give more uniform results, which makes 
  it easier to code for situation where you may use both.

* `phylo_path` and all related methods now deal automatically with both 
  continuous and binary data. All separate binary functions and methods have
  disappeared as they are no longer needed. Mixing of binary and continious
  data in the same models is now allowed.
  
* The variable order in d-seperation statements now better follows the causal
  flow of the DAG.

* Added `plot()` method for `phylopath.summary` objects, that shows the weights
  and p-values for the different models.
  
* `coef_plot()` gained `error_bar`, `order_by`, `from` and `to` arguments. The 
  first allows the user to choose between confidence invervals and standard 
  errors, the second to order the paths by several methods, and the last two
  can be used to select only certain paths.
  
* Plotting methods of causal models now support a manual layout.

* Plotting of fitted DAG's now uses edge width instead of color to indicate, 
  the standardized regression coefficient strength, but this can be reverted 
  using the `type` argument.

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
  All package versions of dependencies are now defined (reported by 
  @ManuelaGonzalez).

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
  
