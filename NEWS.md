phylopath 1.4.0
--------------------------------------------------------------------------------

This release is mainly aimed at avoiding a common confusion about `phylopath`
output, that was long-standing and avoidable. When models would use both
continuous and binary variables, the coefficients would be very confusing to
interpret, since they mix scales and are not consistently standardized. The
package now tries to be much more vocal about what path coefficients actually
mean and warn the user about mixing coefficients with different types.

Note that this concerns the coefficients only, comparison of the causal models
is not an issue.

New features in this release:

* Models that contain both binary and continuous variables now warn, once per
  session, that their coefficients cannot be compared with one another. Every plot
  and printout of such a model repeats this more briefly, and `coef_plot()` refuses
  to sort the paths by `"strength"`, since that ordering would be partly an
  artefact of the scales.

* `plot()` on a fitted model no longer labels its colour legend "standardized path
  coefficient" regardless of what the coefficients are, but names the scale that is
  actually shown.

* `fitted_DAG` objects gain a `binary` element, recording which variables were
  modelled as binary. Objects made by earlier versions do not have it, and are
  reported without a scale rather than being assumed to be standardized.

* New `as.data.frame()` method for fitted models, returning the estimated paths
  as a `data.frame` with one row per path. This function supports a `order_by`
  argument just like `coef_plot()` to order the paths.

There are also a lot of bug fixes and minor improvements in reporting. It is
unlikely any of these have affected real analyses, as they deal with either
imprecise reporting or very unlikely edge cases:

* Bug fix: `coef_plot(order_by = "causal")` derived the causal ordering from the
  positive coefficients only, so the ordering was incorrect in the presence of
  negative paths.

* Bug fix: `CICc` could not be calculated when the number of parameters was
  exactly one less than the number of species, but silently returned `Inf`
  instead of `NA`. Affected models now return `NA`.

* Bug fix: the `phylo_par` column now correctly reports the alpha parameter in
  the case of GLMs.

* Bug fix: model sets in which only some of the models were named kept empty
  names for the others, which made those models impossible to select with
  `best()` and `choice()`. Unnamed models in a partially named set are now
  labelled in the same way as those in a fully unnamed set.

* `phylo_path` now supports more than 26 models.

* Bug fix: the notice pointing to `show_warnings()` was only shown when more
  than one model produced a warning, so a single warning went unmentioned.

* Bug fix: `show_warnings()` and the error for variables with more than two
  levels appended `FALSE` to their messages.

* `best()` and `average()` now return an informative error when `CICc` could not
  be calculated for one or more of the models, rather than failing with an
  uninformative error.

* Informative errors are now returned when a variable is used in the causal
  models but is not a column of `data`, when the `order` argument does not
  contain exactly the variables used in the models, and when two causal models
  are given the same name.

* All errors, warning and messages have been reviewed for clarity, and many
  now have more consistent language, more information, and give tips on how
  to resolve the issue. They also use `rlang` now, for nicer formatting.

* `est_DAG()` now defaults to `model = "lambda"` and `method = "logistic_MPLE"`,
  matching `phylo_path()`. Previously these arguments had no defaults, so calling
  `est_DAG()` without them failed.

* Bug fix: implemented more consistent handling of very small p-values and
  avoidance of numerical underflow.

* Updated `ggplot2` code to remove deprecated functionality.

* Fixed several inaccuracies in the vignettes, including a reference to a
  function that does not exist, and restored two `coef_plot()` examples in the
  introduction that had been commented out.

* Fitted models now report what their coefficients actually are. Paths into a
  binary variable are log odds ratios, while paths into a continuous variable are
  standardized regression coefficients, and previously both were presented as
  though they were the same quantity.

* New `print()` method for fitted models, which lists the estimated paths with
  their coefficients and intervals, and states which variables are continuous and
  which are binary, rather than printing the underlying matrices.

* `coef_plot()` now names the scale of the coefficients on its axis, instead of
  always describing them as standardized.

* Bug fix: `est_DAG()` did not recognize binary variables supplied as character
  vectors, although those are documented as acceptable, and `phylo_path()` accepts
  them. Such a variable failed with an error from `phylolm` when it was the
  outcome of a path, and was silently treated as continuous when it was the
  predictor. Character variables are now converted to factors, as they already
  were in `phylo_path()`.

* Bug fix: `plot_model_set()` located the edges of a causal model incorrectly
  when a model was not topologically sorted, which could silently draw edges
  that were not in the model and omit ones that were.

* Bug fix: passing an argument that neither `phylolm()` nor `phyloglm()` accepts
  produced one warning per fitted regression, and did not say which argument was
  at fault. Such arguments are now reported once per call, by name.

* The vignettes now have descriptive titles in the vignette index, and all
  figures have alt-text.

* Added a unit test suite covering the d-separation basis set, the model
  comparison statistics, model averaging, and the plotting functions.

phylopath 1.3.1
--------------------------------------------------------------------------------

* Maintenance update to confirm to update CRAN documentation guidelines.

phylopath 1.3.0
--------------------------------------------------------------------------------

* MuMIn was removed as a dependency, as it was threatened to be removed from
  CRAN.

phylopath 1.2.1
--------------------------------------------------------------------------------

* A warning is now generated when a user passes a data column with binary data
  as a numeric vector.

* Informative errors are now returned when the `order` argument does not contain
  each variable exactly once.

* Fixed a rare bug in `find_consensus_order`, due to a particular edge case of
  order combinations. In old R versions this would generate a warning about the
  an `if` condition with length > 1, which in newer versions results in an error.
  (Thanks to Laura Alencar for the report.)

phylopath 1.2.0
--------------------------------------------------------------------------------

* Replaced parallel processing based on the `parallel` package and `pbapply` to
  use `future` instead. Use e.g. `future::plan("multisession", workers = n)` to
  enable parallel processing for both model comparison (parallel over dsep
  statements) and model estimation (when using bootstrapping).

* Fixed a bug that no longer allowed parallel processing in `phylo_path`.

* Fixed a bug where the range of the width scale for paths in `plot.fitted_DAG`
  was incorrectly set to the `max(weight)`, instead of `max(abs(weight))`.
  (Thanks Yu Xu for the report.)

* Better user messaging and documentation around the use of the `boot` parameter.

* Update of binary vignette to include more info on convergence warnings.

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
