# The small internal helpers: capturing model output, extracting quantities from
# fitted models, and routing arguments to the right fitting function.

test_that("quiet_safely captures results, errors and warnings separately", {
  ok <- quiet_safely(function() 42)()
  expect_equal(ok$result, 42)
  expect_null(ok$error)
  expect_null(ok$warning)

  failed <- quiet_safely(function() stop("boom"))()
  expect_null(failed$result)
  expect_equal(failed$error, "boom")

  warned <- quiet_safely(function() {
    warning("careful")
    1
  })()
  expect_equal(warned$result, 1)
  expect_equal(warned$warning, "careful")
  expect_null(warned$error)

  # Warnings must not be emitted to the user, they are collected instead.
  expect_no_warning(quiet_safely(function() { warning("quiet"); 1 })())

  # Several warnings are all kept.
  many <- quiet_safely(function() {
    warning("one")
    warning("two")
    3
  })()
  expect_equal(many$warning, c("one", "two"))

  # A warning raised before an error keeps both.
  both <- quiet_safely(function() {
    warning("first")
    stop("then")
  })()
  expect_equal(both$error, "then")
  expect_equal(both$warning, "first")
})

test_that("phylo_g_lm chooses phylolm for numeric and phyloglm for factor responses", {
  data <- test_data()
  tree <- test_tree()

  cont <- phylo_g_lm(B ~ A, data, tree, model = "BM", method = "logistic_MPLE")
  expect_s3_class(cont$result, "phylolm")
  expect_null(cont$error)

  data$bin <- factor(rep(c("no", "yes"), 20))
  bin <- phylo_g_lm(bin ~ A, data, tree, model = "BM", method = "logistic_MPLE")
  expect_s3_class(bin$result, "phyloglm")
})

test_that("phylo_g_lm converts a factor response to 0/1 using the level order", {
  data <- test_data()
  # A binary response perfectly determined by a threshold on A, so the sign of
  # the coefficient tells us which level was coded as 1.
  data$bin <- factor(ifelse(data$A > stats::median(data$A), "high", "low"),
                     levels = c("low", "high"))
  fit <- phylo_g_lm(bin ~ A, data, test_tree(), model = "BM",
                    method = "logistic_MPLE")$result
  # "low" is the first level, so it is coded 0 and "high" is 1; higher A must
  # therefore raise the fitted probability.
  expect_gt(stats::coef(fit)[["A"]], 0)

  # Reversing the level order reverses the sign.
  data$bin <- factor(data$bin, levels = c("high", "low"))
  fit2 <- phylo_g_lm(bin ~ A, data, test_tree(), model = "BM",
                     method = "logistic_MPLE")$result
  expect_lt(stats::coef(fit2)[["A"]], 0)
})

test_that("phylo_g_lm warns about arguments it cannot route", {
  expect_warning(
    phylo_g_lm(B ~ A, test_data(), test_tree(), model = "BM",
               method = "logistic_MPLE", dots = list(nonsense = 1)),
    "not recognized"
  )
  # A genuine phylolm argument is accepted quietly.
  expect_no_warning(
    phylo_g_lm(B ~ A, test_data(), test_tree(), model = "lambda",
               method = "logistic_MPLE", dots = list(lower.bound = 0.01))
  )
})

test_that("phylo_g_lm strips the call from the fitted model", {
  # quiet_safely mangles the call, so it is removed to keep printing readable.
  fit <- phylo_g_lm(B ~ A, test_data(), test_tree(), model = "BM",
                    method = "logistic_MPLE")
  expect_null(fit$result$call)
})

test_that("get_p reads the last coefficient and floors at machine epsilon", {
  fit <- phylolm::phylolm(B ~ A, data = test_data(), phy = test_tree(), model = "BM")
  tbl <- stats::coef(phylolm::summary.phylolm(fit))
  expect_equal(get_p(fit), tbl[nrow(tbl), "p.value"])

  # A p-value of exactly 0 would make Fisher's C infinite, so it is floored.
  fake <- fit
  fake$coefficients[2] <- 1e10
  expect_gte(get_p(fake), .Machine$double.eps)
  expect_true(is.finite(C_stat(get_p(fake))))
})

test_that("get_est and get_se drop the intercept", {
  fit <- phylolm::phylolm(C ~ A + B, data = test_data(), phy = test_tree(), model = "BM")
  tbl <- stats::coef(phylolm::summary.phylolm(fit))
  expect_equal(get_est(fit), stats::coef(fit)[-1])
  expect_length(get_est(fit), 2)
  expect_equal(get_se(fit), tbl[-1, "StdErr"])
})

test_that("get_lower and get_upper are NA without bootstrapping", {
  fit <- phylolm::phylolm(B ~ A, data = test_data(), phy = test_tree(), model = "BM")
  expect_true(is.na(get_lower(fit)))
  expect_true(is.na(get_upper(fit)))

  boot_fit <- phylolm::phylolm(B ~ A, data = test_data(), phy = test_tree(),
                               model = "BM", boot = 20)
  expect_false(is.na(get_lower(boot_fit)))
  expect_lt(get_lower(boot_fit), get_upper(boot_fit))
})

test_that("get_phylo_param returns the optimised parameter, or NA if there is none", {
  lambda <- phylolm::phylolm(B ~ A, data = test_data(), phy = test_tree(),
                             model = "lambda")
  expect_equal(get_phylo_param(lambda), lambda$optpar)
  expect_false(is.na(get_phylo_param(lambda)))

  # Brownian motion has no free parameter to report.
  bm <- phylolm::phylolm(B ~ A, data = test_data(), phy = test_tree(), model = "BM")
  expect_null(bm$optpar)
  expect_true(is.na(get_phylo_param(bm)))
})

test_that("combine_with_labels relabels nodes and validates the input", {
  l <- data.frame(name = c("a", "b"), x = 1:2)
  class(l) <- c("layout_ggraph", "data.frame")
  out <- combine_with_labels(l, c(a = "Alpha", b = "Beta"))
  expect_equal(as.character(out$name), c("Alpha", "Beta"))
  # The layout class is preserved, or ggraph would not know what it has.
  expect_s3_class(out, "layout_ggraph")

  expect_equal(combine_with_labels(l, NULL), l)
  expect_error(combine_with_labels(l, c("Alpha", "Beta")), "must be a named vector")
  expect_error(combine_with_labels(l, c(a = "Alpha")), "missing from labels")
})

test_that("show_warnings requires a phylopath object", {
  expect_error(show_warnings(structure(list(), class = "nope")),
               "expects a phylopath object")
})

test_that("show_warnings re-emits the collected warnings", {
  fake <- structure(list(warnings = list("first problem", "second problem")),
                    class = "phylopath")
  # Both warnings are re-emitted, so capture them in turn.
  expect_warning(expect_warning(show_warnings(fake), "first problem"), "second problem")
  # It returns nothing useful, it is called for its side effect.
  expect_null(suppressWarnings(show_warnings(fake)))
})

test_that("an analysis without warnings has an empty warning list", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  expect_length(p$warnings, 0)
  expect_no_warning(show_warnings(p))
})

test_that("get_phylo_param reports alpha for phyloglm models", {
  # phyloglm stores its phylogenetic parameter as `alpha`; reading only `optpar`
  # made phylo_par NA for every binary variable.
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  fit <- suppressWarnings(phylo_g_lm(bin ~ A, data, test_tree(), model = "BM",
                                     method = "logistic_MPLE")$result)
  expect_s3_class(fit, "phyloglm")
  expect_null(fit$optpar)
  expect_false(is.null(fit$alpha))
  expect_equal(get_phylo_param(fit), fit$alpha)
  expect_false(is.na(get_phylo_param(fit)))
})

test_that("phylo_par is populated for an all-binary analysis", {
  p <- suppressMessages(suppressWarnings(phylo_path(
    define_model_set(a = c(C ~ M, P ~ C), .common = c(G ~ G, D ~ D)),
    cichlids, cichlids_tree
  )))
  expect_false(any(is.na(p$d_sep$a$phylo_par)))
  expect_true(all(p$d_sep$a$phylo_par > 0))
})

