# Model averaging. The expected values here are computed on paper from the
# definitions in the documentation and von Hardenberg & Gonzalez-Voyer (2013):
# full averaging treats a missing path as a zero, conditional averaging averages
# only over the models where the path appears.

# Two models over A, B, C:
#   m1: A -> B = 1.0 (se 0.1), B -> C = 0.5 (se 0.2)
#   m2: A -> B = 2.0 (se 0.3)          [no B -> C]
# with weights 0.75 and 0.25.
avg_fixture <- function() {
  list(
    m1 = make_fitted_DAG(c("A", "B", "C"),
                         list(AB = 1.0, BC = 0.5), list(AB = 0.1, BC = 0.2)),
    m2 = make_fitted_DAG(c("A", "B", "C"),
                         list(AB = 2.0), list(AB = 0.3))
  )
}

test_that("a path present in all models is a plain weighted mean", {
  f <- avg_fixture()
  for (method in c("full", "conditional")) {
    avg <- average_DAGs(f, c(0.75, 0.25), method)
    expect_equal(avg$coef["A", "B"], 0.75 * 1.0 + 0.25 * 2.0)
  }
})

test_that("full averaging shrinks a path that is missing from some models", {
  avg <- average_DAGs(avg_fixture(), c(0.75, 0.25), "full")
  # The missing coefficient counts as zero: 0.75 * 0.5 + 0.25 * 0.
  expect_equal(avg$coef["B", "C"], 0.375)
  # and the SE accounts for both within model variance and the spread.
  expect_equal(
    avg$se["B", "C"],
    sqrt(0.75 * (0.2^2 + (0.5 - 0.375)^2) + 0.25 * (0^2 + (0 - 0.375)^2))
  )
})

test_that("conditional averaging ignores models where the path is absent", {
  avg <- average_DAGs(avg_fixture(), c(0.75, 0.25), "conditional")
  # Only m1 has B -> C, so the average is m1's estimate untouched.
  expect_equal(avg$coef["B", "C"], 0.5)
  expect_equal(avg$se["B", "C"], 0.2)
})

test_that("conditional averaging gives larger estimates than full averaging", {
  # This is the documented trade off: conditional averaging is more sensitive to
  # effects that appear in only some models.
  f <- avg_fixture()
  cond <- average_DAGs(f, c(0.75, 0.25), "conditional")
  full <- average_DAGs(f, c(0.75, 0.25), "full")
  expect_gt(abs(cond$coef["B", "C"]), abs(full$coef["B", "C"]))
})

test_that("paths absent from every model stay exactly zero", {
  for (method in c("full", "conditional")) {
    avg <- average_DAGs(avg_fixture(), c(0.75, 0.25), method)
    expect_equal(avg$coef["A", "C"], 0)
    expect_equal(avg$se["A", "C"], 0)
    expect_equal(avg$lower["A", "C"], 0)
    expect_equal(avg$upper["A", "C"], 0)
  }
})

test_that("confidence intervals are normal intervals around the average", {
  avg <- average_DAGs(avg_fixture(), c(0.75, 0.25), "conditional")
  expect_equal(avg$lower, avg$coef - stats::qnorm(0.975) * avg$se)
  expect_equal(avg$upper, avg$coef + stats::qnorm(0.975) * avg$se)
})

test_that("weights are renormalised, so only their ratios matter", {
  f <- avg_fixture()
  expect_equal(average_DAGs(f, c(3, 1), "conditional"),
               average_DAGs(f, c(0.75, 0.25), "conditional"))
  expect_equal(average_DAGs(f, c(30, 10), "full"),
               average_DAGs(f, c(0.75, 0.25), "full"))
})

test_that("equal weights are the default", {
  f <- avg_fixture()
  expect_equal(average_DAGs(f), average_DAGs(f, c(0.5, 0.5)))
  expect_equal(average_DAGs(f)$coef["A", "B"], 1.5)
})

test_that("averaging a single model returns that model's estimates", {
  f <- avg_fixture()
  avg <- average_DAGs(f[1], 1, "conditional")
  expect_equal(avg$coef, f$m1$coef)
  expect_equal(avg$se, f$m1$se)
})

test_that("models whose variables are in different orders are aligned first", {
  # est_DAG returns matrices ordered by each DAG's own topological sort, so the
  # matrices being averaged often disagree on row/column order. Averaging must
  # align them by name, not by position.
  f <- avg_fixture()
  shuffled <- f$m2
  ord <- c("C", "B", "A")
  shuffled$coef <- shuffled$coef[ord, ord]
  shuffled$se <- shuffled$se[ord, ord]

  avg <- average_DAGs(list(f$m1, shuffled), c(0.75, 0.25), "conditional")
  expect_equal(rownames(avg$coef), rownames(f$m1$coef))
  expect_equal(avg$coef["A", "B"], 1.25)
  expect_equal(avg$coef["B", "C"], 0.5)
})

test_that("avg_method is matched and validated", {
  f <- avg_fixture()
  expect_equal(average_DAGs(f, avg_method = "cond"), average_DAGs(f, avg_method = "conditional"))
  expect_error(average_DAGs(f, avg_method = "nonsense"), "should be one of")
  expect_error(average_DAGs(f, c(1, 1), "full", extra = 1), "deprecated")
})

test_that("the result is a fitted_DAG with all four components", {
  avg <- average_DAGs(avg_fixture())
  expect_s3_class(avg, "fitted_DAG")
  expect_named(avg, c("coef", "se", "lower", "upper", "binary"))
  expect_equal(dimnames(avg$se), dimnames(avg$coef))
})

test_that("averaging real fitted DAGs reproduces a hand-computed weighted mean", {
  data <- test_data()
  tree <- test_tree()
  d1 <- DAG(B ~ A, C ~ B)
  d2 <- DAG(B ~ A, C ~ A)
  f1 <- est_DAG(d1, data, tree, model = "BM", method = "logistic_MPLE")
  f2 <- est_DAG(d2, data, tree, model = "BM", method = "logistic_MPLE")
  avg <- average_DAGs(list(f1, f2), c(0.6, 0.4), "full")

  # A -> B is shared, so it is a straight weighted mean of the two estimates.
  expect_equal(avg$coef["A", "B"],
               0.6 * f1$coef["A", "B"] + 0.4 * f2$coef["A", "B"])
  # B -> C only exists in the first model.
  expect_equal(avg$coef["B", "C"], 0.6 * f1$coef["B", "C"])
})

test_that("averaging carries the binary flags through", {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  f1 <- suppressWarnings(est_DAG(DAG(bin ~ A, B ~ A), data, test_tree(), model = "BM"))
  f2 <- suppressWarnings(est_DAG(DAG(bin ~ A, B ~ bin), data, test_tree(), model = "BM"))

  avg <- average_DAGs(list(f1, f2))
  expect_equal(avg$binary, f1$binary)
  # The coefficient matrices keep their shape: `binary` is not reordered as one.
  expect_equal(dim(avg$coef), dim(f1$coef))
  expect_equal(dimnames(avg$coef), dimnames(f1$coef))
})

test_that("averaging refuses models that disagree about which variables are binary", {
  # Averaging a log odds ratio with a standardized coefficient would produce a
  # number that means nothing.
  data <- test_data()
  bin_data <- data
  bin_data$B <- factor(ifelse(data$B > stats::median(data$B), "hi", "lo"),
                       levels = c("lo", "hi"))
  f_cont <- est_DAG(DAG(B ~ A), data, test_tree(), model = "BM")
  f_bin  <- suppressWarnings(est_DAG(DAG(B ~ A), bin_data, test_tree(), model = "BM"))

  expect_error(average_DAGs(list(f_cont, f_bin)), "disagree about which variables are binary")
  expect_error(average_DAGs(list(f_cont, f_bin)), "\\bB\\b")
})

test_that("averaging models made before binary was recorded warns, but works", {
  # The averaged model cannot be labelled either, so this is said at the point
  # where the information is lost, not only when the result is later used.
  f1 <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  f2 <- est_DAG(DAG(B ~ A, C ~ A), test_data(), test_tree(), model = "BM")
  old <- lapply(list(f1, f2), function(f) { f$binary <- NULL; f })

  expect_warning(average_DAGs(old), "Cannot determine which variables")
  avg <- suppressWarnings(average_DAGs(old))
  # Nothing is invented: the flags stay absent, so nothing is claimed downstream.
  expect_null(avg$binary)
  expect_equal(dim(avg$coef), dim(f1$coef))
})
