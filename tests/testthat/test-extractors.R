# best(), choice() and average() pull models out of a phylopath object and
# estimate them. They must select the right model and carry the original
# analysis settings through to the refit.

pp <- function() {
  phylo_path(
    define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A), c = c(C ~ A, B ~ C),
                     .common = c(D ~ D)),
    test_data(), test_tree(), model = "BM"
  )
}

test_that("best returns the model with the lowest CICc, fitted", {
  p <- pp()
  s <- summary(p)
  b <- best(p)

  expect_s3_class(b, "fitted_DAG")
  expect_equal(b, est_DAG(p$model_set[[s$model[1]]], p$data, p$tree, p$model, p$method))
  # The paths of the best model, and only those, are non-zero.
  expect_equal(b$coef != 0, p$model_set[[s$model[1]]] == 1, ignore_attr = TRUE)
})

test_that("choice selects by name and by index", {
  p <- pp()
  by_name <- choice(p, "b")
  by_index <- choice(p, 2)

  expect_equal(by_name, by_index)
  expect_equal(by_name, est_DAG(p$model_set$b, p$data, p$tree, p$model, p$method))
  expect_equal(by_name$coef != 0, p$model_set$b == 1, ignore_attr = TRUE)
  # Choosing the top model is the same as best().
  expect_equal(choice(p, summary(p)$model[1]), best(p))
})

test_that("choice rejects a model that does not exist", {
  expect_error(choice(pp(), "does_not_exist"))
  expect_error(choice(pp(), 99))
})

test_that("the extractors require a phylopath object", {
  expect_error(best(1), "inherits")
  expect_error(choice(1, "a"), "inherits")
  expect_error(average(1), "inherits")
})

test_that("average averages the models within the cut off", {
  p <- pp()
  s <- summary(p)
  # With this fixture one model is far better than the rest, so the default
  # cut off keeps a single model and averaging is a no-op.
  expect_equal(sum(s$delta_CICc < 2), 1)
  # average() tags its coef matrix with the DAG class, so compare values.
  expect_equal(unclass(average(p)$coef), unclass(best(p)$coef))

  # cut_off = Inf averages over everything.
  avg_all <- average(p, cut_off = Inf)
  fits <- lapply(p$model_set[s$model], est_DAG, p$data, p$tree, p$model, p$method)
  expect_equal(unclass(avg_all$coef), unclass(average_DAGs(fits, s$w, "conditional")$coef))
})

test_that("average passes avg_method through", {
  p <- pp()
  s <- summary(p)
  cond <- average(p, cut_off = Inf, avg_method = "conditional")
  full <- average(p, cut_off = Inf, avg_method = "full")
  expect_false(isTRUE(all.equal(cond$coef, full$coef)))

  fits <- lapply(p$model_set[s$model], est_DAG, p$data, p$tree, p$model, p$method)
  expect_equal(unclass(full$coef), unclass(average_DAGs(fits, s$w, "full")$coef))
})

test_that("an averaged model need not be a DAG", {
  # Averaging models that disagree about direction can produce a matrix with
  # paths both ways, which is expected and documented.
  p <- pp()
  avg <- average(p, cut_off = Inf)
  expect_s3_class(avg$coef, "DAG")
  expect_true(sum(avg$coef != 0) > 0)
})

test_that("settings from the original call are reused by the extractors", {
  # The evolutionary model must not silently revert to the default.
  p_bm <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
                     test_data(), test_tree(), model = "BM")
  p_lambda <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
                         test_data(), test_tree(), model = "lambda")
  expect_false(isTRUE(all.equal(best(p_bm)$coef, best(p_lambda)$coef)))
  expect_equal(best(p_bm), est_DAG(p_bm$model_set$a, p_bm$data, p_bm$tree, "BM",
                                   "logistic_MPLE"))
})

test_that("boot can be supplied to the extractors and yields intervals", {
  p <- pp()
  b <- suppressWarnings(best(p, boot = 20))
  expect_named(b, c("coef", "se", "lower", "upper"))
  expect_lt(b$lower["A", "B"], b$coef["A", "B"])
  expect_gt(b$upper["A", "B"], b$coef["A", "B"])
})

test_that("combine_dots lets new arguments win over stored ones", {
  expect_equal(combine_dots(list(a = 1, b = 2)), list(a = 1, b = 2))
  expect_equal(combine_dots(list(a = 1, b = 2), a = 99), list(a = 99, b = 2))
  expect_equal(combine_dots(list(a = 1), c = 3), list(c = 3, a = 1))
  expect_length(combine_dots(list()), 0)
})

test_that("best and average refuse to rank models when CICc is missing", {
  # Previously this surfaced as `inherits(DAG, "DAG") is not TRUE`, from
  # selecting NA rows out of the summary table.
  ms <- define_model_set(sparse = c(B ~ A, C ~ B), dense = c(B ~ A, C ~ B, D ~ C))
  p <- suppressWarnings(phylo_path(ms, small_data(), small_tree(), model = "BM"))

  expect_error(suppressWarnings(best(p)), "CICc could not be calculated")
  expect_error(suppressWarnings(average(p)), "CICc could not be calculated")
  # The message names the model at fault and says what to do about it.
  expect_error(suppressWarnings(best(p)), "dense")
  expect_error(suppressWarnings(best(p)), "n > q \\+ 1")

  # choice() does not rank models, so it still works.
  expect_s3_class(suppressWarnings(choice(p, "sparse")), "fitted_DAG")
})

test_that("the averaged coefficient matrix keeps the modern matrix class", {
  avg <- average(pp(), cut_off = Inf)
  expect_equal(class(avg$coef), c("matrix", "array", "DAG"))
  expect_true(is.matrix(avg$coef))
  expect_true(is.array(avg$coef))
})
