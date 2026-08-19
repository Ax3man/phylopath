# est_DAG() fits the paths of a single DAG and returns standardized path
# coefficients. The published method reports standardized coefficients, so the
# scaling of the data and the placement of coefficients into the matrix both
# matter.

test_that("est_DAG returns matrices shaped like the DAG", {
  d <- DAG(B ~ A, C ~ B)
  fit <- est_DAG(d, test_data(), test_tree(), model = "BM", method = "logistic_MPLE")

  expect_s3_class(fit, "fitted_DAG")
  expect_named(fit, c("coef", "se", "binary"))
  expect_equal(dim(fit$coef), dim(d))
  expect_equal(dimnames(fit$coef), dimnames(d))
  expect_equal(dimnames(fit$se), dimnames(d))
})

test_that("coefficients are standardized, matching phylolm on scaled data", {
  data <- test_data()
  tree <- test_tree()
  fit <- est_DAG(DAG(B ~ A, C ~ B), data, tree, model = "BM", method = "logistic_MPLE")

  scaled <- data
  scaled[] <- lapply(scaled, function(x) as.numeric(scale(x)))
  rownames(scaled) <- rownames(data)

  m_ab <- phylolm::phylolm(B ~ A, data = scaled, phy = tree, model = "BM")
  m_bc <- phylolm::phylolm(C ~ B, data = scaled, phy = tree, model = "BM")

  expect_equal(fit$coef["A", "B"], stats::coef(m_ab)[[2]])
  expect_equal(fit$coef["B", "C"], stats::coef(m_bc)[[2]])
  expect_equal(fit$se["A", "B"],
               stats::coef(phylolm::summary.phylolm(m_ab))[2, "StdErr"])
  expect_equal(fit$se["B", "C"],
               stats::coef(phylolm::summary.phylolm(m_bc))[2, "StdErr"])
})

test_that("scaling makes the coefficients invariant to units", {
  # A standardized coefficient must not change if a variable is rescaled.
  data <- test_data()
  rescaled <- data
  rescaled$A <- rescaled$A * 1000 + 50
  fit1 <- est_DAG(DAG(B ~ A, C ~ B), data, test_tree(), model = "BM",
                  method = "logistic_MPLE")
  fit2 <- est_DAG(DAG(B ~ A, C ~ B), rescaled, test_tree(), model = "BM",
                  method = "logistic_MPLE")
  expect_equal(fit1$coef, fit2$coef)
  expect_equal(fit1$se, fit2$se)
})

test_that("coefficients land in the cell for their own path", {
  # A multi-parent node is the case where a mis-ordered assignment would show
  # up: each parent's estimate must go into its own row.
  data <- test_data()
  tree <- test_tree()
  d <- DAG(C ~ A + B, B ~ A)
  fit <- est_DAG(d, data, tree, model = "BM", method = "logistic_MPLE")

  scaled <- data
  scaled[] <- lapply(scaled, function(x) as.numeric(scale(x)))
  rownames(scaled) <- rownames(data)
  m <- phylolm::phylolm(C ~ A + B, data = scaled, phy = tree, model = "BM")

  expect_equal(fit$coef["A", "C"], stats::coef(m)[["A"]])
  expect_equal(fit$coef["B", "C"], stats::coef(m)[["B"]])
  # Absent paths stay exactly zero.
  expect_equal(fit$coef["C", "A"], 0)
  expect_equal(fit$coef["C", "B"], 0)
  expect_equal(sum(fit$coef != 0), 3)
})

test_that("a node with no parents contributes only zeros", {
  fit <- est_DAG(DAG(B ~ A, C ~ C), test_data(), test_tree(), model = "BM",
                 method = "logistic_MPLE")
  expect_true(all(fit$coef[, "A"] == 0))
  expect_true(all(fit$coef[, "C"] == 0))
  expect_equal(sum(fit$coef != 0), 1)
})

test_that("boot adds bootstrap confidence intervals that bracket the estimate", {
  fit <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                  model = "BM", method = "logistic_MPLE", boot = 30))
  expect_named(fit, c("coef", "se", "lower", "upper", "binary"))
  for (path in list(c("A", "B"), c("B", "C"))) {
    expect_lt(fit$lower[path[1], path[2]], fit$coef[path[1], path[2]])
    expect_gt(fit$upper[path[1], path[2]], fit$coef[path[1], path[2]])
  }
  # Absent paths have no interval.
  expect_equal(fit$lower["B", "A"], 0)
  expect_equal(fit$upper["B", "A"], 0)
})

test_that("est_DAG requires a DAG", {
  expect_error(
    est_DAG(matrix(0, 2, 2), test_data(), test_tree(), model = "BM"),
    "`DAG` must be a DAG object"
  )
  # And says where one comes from, rather than quoting the failed check.
  expect_error(
    est_DAG(matrix(0, 2, 2), test_data(), test_tree(), model = "BM"),
    "define_model_set"
  )
})

test_that("errors in a path model are reported with the offending formula", {
  expect_error(
    est_DAG(DAG(B ~ A), test_data(), test_tree(), model = "not_a_model"),
    "A phylogenetic regression could not be fitted"
  )
})

test_that("binary responses are fitted with phyloglm", {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  fit <- suppressWarnings(est_DAG(DAG(bin ~ A), data, test_tree(), model = "BM",
                                 method = "logistic_MPLE"))
  # A coefficient was estimated for the path into the binary variable.
  expect_false(fit$coef["A", "bin"] == 0)
  expect_false(fit$se["A", "bin"] == 0)
})

test_that("est_DAG has usable defaults for model and method", {
  # Both are documented via @inheritParams, but had no defaults here, so a
  # perfectly reasonable call failed with "argument \"model\" is missing".
  expect_equal(formals(est_DAG)$model, "lambda")
  expect_equal(formals(est_DAG)$method, "logistic_MPLE")

  fit <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree())
  expect_s3_class(fit, "fitted_DAG")
  # The default really is Pagel's lambda, matching phylo_path().
  expect_equal(fit$coef, est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                 model = "lambda")$coef)

  # The documented example, which previously worked only by accident.
  expect_s3_class(est_DAG(DAG(LS ~ BM, NL ~ BM, DD ~ NL + LS), rhino, rhino_tree,
                          "lambda"),
                  "fitted_DAG")
})

test_that("est_DAG defaults work for binary variables too", {
  # method is only consulted for binary responses, which is why the missing
  # default went unnoticed for continuous data.
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  expect_s3_class(suppressWarnings(est_DAG(DAG(bin ~ A), data, test_tree())),
                  "fitted_DAG")
})

test_that("est_DAG records which variables were modelled as binary", {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  fit <- suppressWarnings(est_DAG(DAG(bin ~ A, B ~ A), data, test_tree(), model = "BM"))

  # One entry per variable in the DAG, named, and matching how phylo_g_lm()
  # decides between phylolm and phyloglm.
  expect_named(fit$binary, rownames(fit$coef))
  expect_true(fit$binary[["bin"]])
  expect_false(any(fit$binary[c("A", "B")]))

  # Continuous-only models say so, rather than leaving it unknown.
  cont <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  expect_false(any(cont$binary))
})

test_that("character columns are treated as binary, like phylo_path does", {
  # `data` is documented as accepting binary variables "either as character values
  # or factors", but only factors were recognized: phylo_path() converts them via
  # check_models_data_tree(), while est_DAG() did not. A character response then
  # failed inside phylolm, and a character predictor was silently labelled
  # continuous even though it was fitted as a contrast between its two levels.
  data <- test_data()
  data$chr <- rep(c("no", "yes"), 20)
  data$fac <- factor(data$chr)

  # As a response: fits, rather than failing inside phylolm.
  resp <- suppressWarnings(est_DAG(DAG(chr ~ A), data, test_tree(), model = "BM"))
  expect_true(resp$binary[["chr"]])

  # As a predictor: identical to the factor version, and labelled the same way.
  chr <- suppressWarnings(est_DAG(DAG(B ~ chr), data, test_tree(), model = "BM"))
  fac <- suppressWarnings(est_DAG(DAG(B ~ fac), data, test_tree(), model = "BM"))
  expect_true(chr$binary[["chr"]])
  expect_equal(unname(chr$coef["chr", "B"]), unname(fac$coef["fac", "B"]))
  # Which means such a model is now recognized as mixed, and so gets labelled.
  expect_true(is_mixed(chr))
})

test_that("est_DAG rejects variables that are not binary", {
  data <- test_data()
  data$three <- rep(c("a", "b", "c"), length.out = 40)
  expect_error(est_DAG(DAG(B ~ three), data, test_tree(), model = "BM"),
               "`three` is expected to be binary, but has 3 levels")
})

test_that("est_DAG ignores columns that are not in the DAG", {
  # `data` may hold other columns, which are none of est_DAG's business, so a
  # many-levelled column elsewhere in the data must not be an error.
  data <- test_data()
  data$unused <- rep(c("a", "b", "c"), length.out = 40)
  expect_no_error(est_DAG(DAG(B ~ A), data, test_tree(), model = "BM"))
})
