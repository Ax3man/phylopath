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

# A fitted_DAG whose responses are of both types, so its coefficients are on two
# incomparable scales.
mixed_scales <- function() {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  suppressWarnings(est_DAG(DAG(bin ~ A, C ~ bin), data, test_tree(), model = "BM"))
}

# A fitted_DAG whose every variable is binary, so all its coefficients are level to
# level log odds ratios and are comparable with each other.
all_binary <- function() {
  data <- test_data()
  data$bin1 <- factor(rep(c("no", "yes"), 20))
  data$bin2 <- factor(rep(c("lo", "hi"), each = 20))
  suppressWarnings(est_DAG(DAG(bin2 ~ bin1), data, test_tree(), model = "BM"))
}

test_that("coef_plot names the scale on the axis when there is only one", {
  cont <- coef_plot(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM"),
                    error_bar = "se")
  expect_equal(cont$labels$y, "standardized regression coefficient \U00B1 SE")

  # An all-binary model: both ends of the path are binary, so every coefficient is
  # a level to level log odds ratio.
  expect_equal(coef_plot(all_binary(), error_bar = "se")$labels$y,
               "log odds ratio \U00B1 SE")
})

test_that("coef_plot claims no scale for a mixed model", {
  # The paths stay on one panel, but the axis no longer asserts a scale that only
  # some of the points are on. The caption and the warning cover the rest.
  p <- suppressWarnings(coef_plot(mixed_scales(), error_bar = "se"))
  expect_equal(p$labels$y, "path coefficient \U00B1 SE")
  expect_no_error(ggplot2::ggplot_build(p))
  # The vignettes flip the coordinates, so that combination has to survive.
  expect_no_error(ggplot2::ggplot_build(
    suppressWarnings(coef_plot(mixed_scales(), error_bar = "se", reverse_order = TRUE)) +
      ggplot2::coord_flip()))
})

test_that("an unknown coefficient scale warns wherever it is reported", {
  # Objects made before the binary flags existed get the neutral label rather
  # than the old one, which asserted "standardized" whether or not that was true,
  # and every function that reports coefficients says why it cannot label them.
  fit <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  fit$binary <- NULL

  expect_warning(coef_plot(fit, error_bar = "se"), "Cannot determine which variables")
  expect_warning(plot(fit), "Cannot determine which variables")
  expect_warning(capture.output(print(fit)), "Cannot determine which variables")
  # The warning explains what the two possible scales are, and how to fix it.
  expect_warning(coef_plot(fit, error_bar = "se"), "log odds ratios")
  expect_warning(coef_plot(fit, error_bar = "se"), "Refit the model")

  p <- suppressWarnings(coef_plot(fit, error_bar = "se"))
  expect_equal(p$labels$y, "path coefficient \U00B1 SE")
})

test_that("a mixed model warns once per session, wherever it is reported", {
  # The warning is deliberately once per session: a single analysis fits, prints
  # and plots the same model, and repeating it at each step would be noise.
  f <- mixed_scales()

  rlang::reset_warning_verbosity("phylopath_mixed_type")
  expect_warning(plot(f), "cannot be compared with one another")
  # Already said once, so these stay quiet.
  expect_no_warning(plot(f))
  expect_no_warning(coef_plot(f, error_bar = "se"))
  expect_no_warning(capture.output(print(f)))

  # Each entry point raises it when it is the first to be reached.
  for (fun in list(function() coef_plot(f, error_bar = "se"),
                   function() capture.output(print(f)),
                   function() plot(f))) {
    rlang::reset_warning_verbosity("phylopath_mixed_type")
    expect_warning(fun(), "cannot be compared with one another")
  }

  # A single-scale model never warns.
  rlang::reset_warning_verbosity("phylopath_mixed_type")
  expect_no_warning(plot(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                 model = "BM")))
})

test_that("mixed models are annotated on every output that shows coefficients", {
  # The polite restatement, which is repeated even after the warning has been
  # said, because a plot or printout may be read on its own.
  f <- mixed_scales()
  expect_match(suppressWarnings(coef_plot(f, error_bar = "se"))$labels$caption,
               "cannot be compared")
  expect_match(suppressWarnings(plot(f))$labels$caption, "cannot be compared")
  expect_match(suppressWarnings(capture.output(print(f))), "cannot be compared",
               all = FALSE)

  # And absent when there is nothing to warn about.
  cont <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  expect_null(coef_plot(cont, error_bar = "se")$labels$caption)
  expect_null(plot(cont)$labels$caption)
})

test_that("plot.fitted_DAG names the scale in its colour legend", {
  # The legend title was hard coded to "standardized path coefficient", which is
  # wrong whenever the response is binary.
  legend_title <- function(p) {
    edge <- vapply(p$scales$scales, function(s) 'edge_colour' %in% s$aesthetics,
                   logical(1))
    p$scales$scales[[which(edge)]]$name
  }
  expect_match(legend_title(plot(all_binary(), type = "color")), "log odds ratio")

  cont <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  expect_match(legend_title(plot(cont, type = "color")), "standardized")
})

test_that("coef_plot refuses to sort a mixed model by strength", {
  # Binary predictors are inflated by at least a factor of two, so sorting by
  # size would order the paths partly by arithmetic rather than by biology.
  expect_error(coef_plot(mixed_scales(), error_bar = "se", order_by = "strength"),
               "not meaningful for a model that contains both binary and continuous")
  # The other orderings remain available, and are suggested by the error.
  expect_no_error(suppressWarnings(
    coef_plot(mixed_scales(), error_bar = "se", order_by = "causal")))
  # Sorting by strength is fine when the scales are comparable.
  expect_no_error(coef_plot(est_DAG(DAG(C ~ A + B, B ~ A), test_data(), test_tree(),
                                    model = "BM"),
                            error_bar = "se", order_by = "strength"))
})

test_that("plot_model_set facets over the models", {
  ms <- define_model_set(one = c(B ~ A, C ~ B), two = c(B ~ A, C ~ A))
  p <- plot_model_set(ms)
  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(p$data$name), c("A", "B", "C"))
  expect_s3_class(plot_model_set(ms, algorithm = "sugiyama"), "ggplot")
  expect_s3_class(plot_model_set(ms, nrow = 2), "ggplot")
})

# The arrows plot_model_set() decided to draw, as a "from->to" character vector,
# per model. The layout keeps the graph it was built from as an attribute.
drawn_edges <- function(p) {
  g <- attr(p$data, "graph")
  e <- igraph::as_edgelist(g)
  split(paste0(e[, 1], "->", e[, 2]), igraph::E(g)$model)
}

test_that("plot_model_set draws exactly the edges the models contain", {
  # An unsorted model can have an edge in the last row of its matrix, which is
  # where locating edges by index arithmetic goes wrong.
  m <- DAG(A ~ B, B ~ C, order = FALSE)          # the chain C -> B -> A
  expect_setequal(drawn_edges(plot_model_set(list(one = m)))$one, c("B->A", "C->B"))
  expect_setequal(drawn_edges(plot_model_set(list(one = DAG(A ~ B, order = FALSE))))$one,
                  "B->A")

  ms <- define_model_set(one = c(B ~ A, C ~ B), two = c(B ~ A, C ~ A))
  drawn <- drawn_edges(plot_model_set(ms))
  expect_setequal(drawn$one, c("A->B", "B->C"))
  expect_setequal(drawn$two, c("A->B", "A->C"))
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

test_that("every plot builds, not only constructs", {
  # Variables referenced in aes() are resolved when a plot is *built*, not when
  # it is constructed, so a mis-named column would pass every test above. These
  # build each plot for real.
  f <- fitted()
  fb <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                 model = "BM", method = "logistic_MPLE", boot = 20))
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  plots <- list(
    dag        = plot(DAG(B ~ A, C ~ B)),
    dag_labels = plot(DAG(B ~ A, C ~ B), labels = c(A = "Ay", B = "Bee", C = "Cee")),
    fitted_w   = plot(f, type = "width"),
    fitted_c   = plot(f, type = "color"),
    coef_se    = coef_plot(f, error_bar = "se"),
    coef_ci    = coef_plot(fb),
    coef_str   = coef_plot(f, error_bar = "se", order_by = "strength"),
    model_set  = plot_model_set(define_model_set(one = c(B ~ A, C ~ B),
                                                two = c(B ~ A, C ~ A))),
    summary    = plot(summary(p))
  )
  for (nm in names(plots)) {
    expect_no_error(ggplot2::ggplot_build(plots[[nm]]))
  }
})
