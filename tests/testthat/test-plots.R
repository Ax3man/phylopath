# Plotting. These are deliberately structural rather than visual: they check
# that a ggplot object is produced, that the data going into it is right, and
# that the documented arguments and error conditions behave. Comparing rendered
# output would need vdiffr and would be fragile across ggplot2 versions.

fitted <- function() {
  est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM",
          method = "logistic_MPLE")
}

test_that("plot.DAG produces a ggplot", {
  expect_s3_class(plot(DAG(A ~ B + C)), "ggplot")
  expect_s3_class(plot(DAG(B ~ A, C ~ B), algorithm = "kk"), "ggplot")
})

test_that("plot.DAG accepts a manual layout and node labels", {
  d <- DAG(a ~ b + c + d)
  ml <- data.frame(name = c("a", "b", "c", "d"), x = c(1, 1, 2, 2), y = c(1, 2, 1, 2))
  p <- plot(d, manual_layout = ml)
  expect_s3_class(p, "ggplot")
  # The supplied coordinates are the ones used.
  expect_setequal(p$data$x, ml$x)
  expect_setequal(p$data$y, ml$y)

  labs <- c(a = "Alpha", b = "Beta", c = "Gamma", d = "Delta")
  expect_setequal(as.character(plot(d, labels = labs)$data$name), unname(labs))
})

test_that("labels must be a complete named vector", {
  d <- DAG(a ~ b)
  expect_error(plot(d, labels = c("Alpha", "Beta")), "labels must be a named vector")
  expect_error(plot(d, labels = c(a = "Alpha")), "Some nodes are missing from labels")
})

test_that("a graph with no paths cannot be plotted", {
  expect_error(plot(DAG(A ~ A)), "no paths to plot")
  empty <- fitted()
  empty$coef[] <- 0
  expect_error(plot(empty), "no paths to plot")
})

test_that("rotation and flipping transform the layout", {
  d <- DAG(a ~ b + c + d)
  ml <- data.frame(name = c("a", "b", "c", "d"), x = c(1, 1, 2, 2), y = c(1, 2, 1, 2))
  base <- plot(d, manual_layout = ml)
  expect_equal(plot(d, manual_layout = ml, flip_x = TRUE)$data$x, -base$data$x)
  expect_equal(plot(d, manual_layout = ml, flip_y = TRUE)$data$y, -base$data$y)
  # A 180 degree rotation negates both coordinates.
  rot <- plot(d, manual_layout = ml, rotation = 180)
  expect_equal(rot$data$x, -base$data$x, tolerance = 1e-8)
  expect_equal(rot$data$y, -base$data$y, tolerance = 1e-8)
})

test_that("adjust_layout rotates by the requested angle", {
  l <- data.frame(x = c(1, 0), y = c(0, 1))
  # Rotation is clockwise: (1, 0) goes to (0, -1) at 90 degrees.
  r90 <- adjust_layout(l, 90, FALSE, FALSE)
  expect_equal(r90$x, c(0, 1), tolerance = 1e-8)
  expect_equal(r90$y, c(-1, 0), tolerance = 1e-8)
  # Properties that hold whatever the handedness: a full turn is the identity,
  # a half turn negates both coordinates, and distances are preserved.
  expect_equal(adjust_layout(l, 360, FALSE, FALSE)$x, l$x, tolerance = 1e-8)
  expect_equal(adjust_layout(l, 360, FALSE, FALSE)$y, l$y, tolerance = 1e-8)
  expect_equal(adjust_layout(l, 180, FALSE, FALSE)$x, -l$x, tolerance = 1e-8)
  for (angle in c(37, 90, 145)) {
    rot <- adjust_layout(l, angle, FALSE, FALSE)
    expect_equal(sqrt(rot$x^2 + rot$y^2), sqrt(l$x^2 + l$y^2), tolerance = 1e-8)
  }
})

test_that("plot.fitted_DAG works for both weight types", {
  f <- fitted()
  expect_s3_class(plot(f), "ggplot")
  expect_s3_class(plot(f, type = "color"), "ggplot")
  expect_s3_class(plot(f, type = "colour"), "ggplot")
  expect_error(plot(f, type = "nonsense"), "should be one of")
  expect_warning(plot(f, width_const = 2), "deprecated")
})

test_that("coef_plot shows one row per estimated path", {
  f <- fitted()
  p <- coef_plot(f, error_bar = "se")
  expect_s3_class(p, "ggplot")
  # Two paths in the DAG, so two points.
  expect_equal(nrow(p$data), 2)
  expect_setequal(as.character(p$data$path), c("A \U2192 B", "B \U2192 C"))
  # Error bars are the standard errors.
  expect_equal(p$data$lower, p$data$coef - c(f$se["A", "B"], f$se["B", "C"]))
  expect_equal(p$data$upper, p$data$coef + c(f$se["A", "B"], f$se["B", "C"]))
})

test_that("coef_plot falls back to standard errors when there are no intervals", {
  expect_message(coef_plot(fitted()), "does not contain confidence intervals")
})

test_that("coef_plot uses bootstrap intervals when available", {
  f <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                model = "BM", method = "logistic_MPLE", boot = 20))
  p <- expect_no_message(coef_plot(f))
  expect_equal(p$data$lower, c(f$lower["A", "B"], f$lower["B", "C"]))
  expect_equal(p$data$upper, c(f$upper["A", "B"], f$upper["B", "C"]))
})

test_that("coef_plot can filter and reorder paths", {
  f <- est_DAG(DAG(C ~ A + B, B ~ A), test_data(), test_tree(), model = "BM",
               method = "logistic_MPLE")
  expect_equal(nrow(coef_plot(f, error_bar = "se", from = "A")$data), 2)
  expect_equal(nrow(coef_plot(f, error_bar = "se", to = "C")$data), 2)
  expect_equal(nrow(coef_plot(f, error_bar = "se", from = "A", to = "C")$data), 1)

  by_strength <- coef_plot(f, error_bar = "se", order_by = "strength")$data
  expect_false(is.unsorted(by_strength$coef))

  # Reversing flips the factor levels used for the axis.
  fwd <- coef_plot(f, error_bar = "se")$data
  rev <- coef_plot(f, error_bar = "se", reverse_order = TRUE)$data
  expect_equal(levels(rev$path), rev(levels(fwd$path)))

  expect_error(coef_plot(f, order_by = "nonsense"), "should be one of")
  expect_error(coef_plot(f, error_bar = "nonsense"), "should be one of")
  expect_error(coef_plot(DAG(B ~ A)), "inherits")
})

test_that("plot_model_set facets over the models", {
  ms <- define_model_set(one = c(B ~ A, C ~ B), two = c(B ~ A, C ~ A))
  p <- plot_model_set(ms)
  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(p$data$name), c("A", "B", "C"))
  expect_s3_class(plot_model_set(ms, algorithm = "sugiyama"), "ggplot")
  expect_s3_class(plot_model_set(ms, nrow = 2), "ggplot")
})

test_that("plot_model_set names unnamed model sets and validates input", {
  expect_s3_class(plot_model_set(list(DAG(B ~ A), DAG(A ~ B))), "ggplot")
  expect_error(plot_model_set(list(DAG(B ~ A), "not a DAG")),
               "model_set should be a list of DAG objects")
  expect_error(plot_model_set(DAG(B ~ A)), "list of DAG objects")
  expect_error(plot_model_set(list(DAG(B ~ A), DAG(C ~ A))),
               "All causal models need to include the same variables")
})

test_that("plot.phylopath_summary produces a ggplot", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  expect_s3_class(plot(summary(p)), "ggplot")
  expect_s3_class(plot(summary(p), cut_off = 5), "ggplot")
})
