# print.fitted_DAG() lists the realised paths and their coefficients, states which
# variables are continuous and which are binary, and names the scale of the
# coefficients whenever every one of them is on the same scale.

# A fitted_DAG with a binary response (bin) that is also used as a predictor.
mixed_fit <- function(boot = 0) {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  suppressWarnings(est_DAG(DAG(bin ~ A, C ~ bin), data, test_tree(), model = "BM",
                           boot = boot))
}

test_that("print.fitted_DAG lists the paths and returns its input invisibly", {
  out <- suppressWarnings(capture.output(res <- withVisible(print(mixed_fit()))))
  expect_false(res$visible)
  expect_s3_class(res$value, "fitted_DAG")

  expect_match(out, "2 paths", all = FALSE)
  # All paths in one table, whatever their scale.
  expect_match(out, "A \U2192 bin", all = FALSE)
  expect_match(out, "bin \U2192 C", all = FALSE)
  expect_length(grep("\U2192", out), 2)
})

test_that("print.fitted_DAG names the scale when every coefficient shares one", {
  cont <- capture.output(print(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                      model = "BM")))
  expect_match(cont, "Paths \U2014 standardized regression coefficients", all = FALSE)
  expect_false(any(grepl("log odds", cont)))

  data <- test_data()
  data$bin1 <- factor(rep(c("no", "yes"), 20))
  data$bin2 <- factor(rep(c("lo", "hi"), each = 20))
  bin <- capture.output(print(suppressWarnings(
    est_DAG(DAG(bin2 ~ bin1), data, test_tree(), model = "BM"))))
  expect_match(bin, "Paths \U2014 log odds ratios", all = FALSE)
  expect_false(any(grepl("standardized", bin)))
})

test_that("print.fitted_DAG names no scale for a mixed model", {
  out <- suppressWarnings(capture.output(print(mixed_fit())))
  # A bare heading, since no single scale covers every path.
  expect_match(out, "^Paths\\s*$", all = FALSE)
  expect_false(any(grepl("Paths \U2014", out)))
})

test_that("print.fitted_DAG reports which variables are of which type", {
  # Using the same vocabulary as print.phylopath(), so a reader can work out what
  # any individual path means.
  mixed <- suppressWarnings(capture.output(print(mixed_fit())))
  expect_match(mixed, "Continuous:\\s+A C", all = FALSE)
  expect_match(mixed, "Binary:\\s+bin", all = FALSE)

  # Categories with no variables in them are left out entirely.
  cont <- capture.output(print(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                      model = "BM")))
  expect_match(cont, "Continuous:", all = FALSE)
  expect_false(any(grepl("Binary:", cont)))
})

test_that("print.fitted_DAG lists only realised paths", {
  # The matrix is mostly zeros; those are absent paths, not estimates of zero.
  fit <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  out <- capture.output(print(fit))
  paths <- grep("\U2192", out, value = TRUE)
  expect_length(paths, 2)
  expect_match(out, "2 paths", all = FALSE)
  # D is in the data but has no paths, so it is never listed.
  expect_false(any(grepl("D", paths)))
})

test_that("print.fitted_DAG shows intervals when they exist and errors otherwise", {
  expect_match(capture.output(print(mixed_fit(boot = 20))), "95% CI", all = FALSE)
  expect_match(capture.output(print(mixed_fit())), "se", all = FALSE)
})

test_that("print.fitted_DAG notes a mixed model, and stays quiet otherwise", {
  expect_match(suppressWarnings(capture.output(print(mixed_fit()))),
               "mixes binary and continuous variables", all = FALSE)

  # An all-binary model is not mixed: every coefficient is a level to level log
  # odds ratio, so they are comparable and there is nothing to note.
  data <- test_data()
  data$bin1 <- factor(rep(c("no", "yes"), 20))
  data$bin2 <- factor(rep(c("lo", "hi"), each = 20))
  all_bin <- suppressWarnings(est_DAG(DAG(bin2 ~ bin1), data, test_tree(), model = "BM"))
  out <- capture.output(print(all_bin))
  expect_false(any(grepl("mixes binary and continuous", out)))
  expect_match(out, "log odds ratios", all = FALSE)
})

test_that("print.fitted_DAG warns and makes no claim when the scale is unknown", {
  fit <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  fit$binary <- NULL
  expect_warning(print(fit), "Cannot determine which variables of this model are binary")
  out <- capture.output(suppressWarnings(print(fit)))
  # Neither scale is asserted, since neither can be known, and the paths are
  # still listed so the coefficients remain readable.
  expect_false(any(grepl("standardized regression coefficients|log odds ratios", out)))
  expect_match(out, "^Paths\\s*$", all = FALSE)
  expect_length(grep("\U2192", out), 2)
})

test_that("print.fitted_DAG copes with a model with no paths at all", {
  fit <- est_DAG(DAG(B ~ A), test_data(), test_tree(), model = "BM")
  fit$coef[] <- 0
  out <- capture.output(print(fit))
  expect_match(out, "0 paths", all = FALSE)
})
