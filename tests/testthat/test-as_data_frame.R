# as.data.frame() turns a fitted_DAG into one row per estimated path. It is both
# a user facing convenience and the single place where the coefficient matrices
# are converted, so print.fitted_DAG() and coef_plot() go through it too.

test_that("as.data.frame returns one row per estimated path", {
  fit <- est_DAG(DAG(C ~ A + B, B ~ A), test_data(), test_tree(), model = "BM")
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_named(df, c("from", "to", "coef", "se"))
  expect_equal(nrow(df), 3)
  expect_setequal(paste(df$from, df$to), c("A B", "A C", "B C"))
  # The values are the ones in the matrices, looked up by name.
  for (i in seq_len(nrow(df))) {
    expect_equal(df$coef[i], fit$coef[df$from[i], df$to[i]])
    expect_equal(df$se[i], fit$se[df$from[i], df$to[i]])
  }
  # Absent paths are omitted rather than reported as zero.
  expect_false(any(df$coef == 0))
  expect_equal(nrow(df), sum(fit$coef != 0))
})

test_that("as.data.frame includes bootstrap intervals when they exist", {
  fit <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                  model = "BM", boot = 20))
  df <- as.data.frame(fit)
  expect_named(df, c("from", "to", "coef", "se", "lower", "upper"))
  expect_true(all(df$lower < df$coef))
  expect_true(all(df$upper > df$coef))
})

test_that("as.data.frame orders the paths as asked", {
  fit <- est_DAG(DAG(C ~ A + B, B ~ A), test_data(), test_tree(), model = "BM")

  # By position in the model, which is the order of the matrix.
  v <- rownames(fit$coef)
  d <- as.data.frame(fit, order_by = "default")
  expect_false(is.unsorted(order(match(d$from, v), match(d$to, v))))

  # By size of the coefficient.
  expect_false(is.unsorted(as.data.frame(fit, order_by = "strength")$coef))

  # All orderings hold the same paths, only rearranged.
  paths <- function(d) paste(d$from, d$to)
  expect_setequal(paths(as.data.frame(fit, order_by = "causal")),
                  paths(as.data.frame(fit, order_by = "default")))

  expect_error(as.data.frame(fit, order_by = "nonsense"), "should be one of")
  # Row names do not carry over from the previous ordering.
  expect_equal(rownames(as.data.frame(fit, order_by = "strength")),
               as.character(seq_len(3)))
})

test_that("causal ordering accounts for negative paths", {
  # The ordering used to be derived from `coef > 0`, so a negative path counted
  # as no path at all and the topological sort saw a different causal model.
  data <- test_data()
  data$B <- -data$B
  fit <- est_DAG(DAG(B ~ A, C ~ B), data, test_tree(), model = "BM")
  expect_lt(fit$coef["A", "B"], 0)

  d <- as.data.frame(fit, order_by = "causal")
  # A is upstream of B, which is upstream of C, whatever the signs.
  expect_equal(paste(d$from, d$to), c("A B", "B C"))
})

test_that("as.data.frame refuses to order a mixed model by strength", {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 20))
  fit <- suppressWarnings(est_DAG(DAG(bin ~ A, C ~ bin), data, test_tree(), model = "BM"))

  expect_error(as.data.frame(fit, order_by = "strength"),
               "not meaningful for a causal model that contains both binary and continuous")
  # The other orderings are still available.
  expect_no_error(as.data.frame(fit, order_by = "causal"))
  expect_no_error(as.data.frame(fit))
})

test_that("as.data.frame copes with a model with no paths", {
  fit <- est_DAG(DAG(B ~ A), test_data(), test_tree(), model = "BM")
  fit$coef[] <- 0
  df <- as.data.frame(fit)
  # An empty table, not NULL, so that callers need no special case.
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0)
  expect_named(df, c("from", "to", "coef", "se"))
})

test_that("causal ordering is correct for an averaged, unsorted model", {
  # Any model whose matrix is not already topologically sorted will order
  # differently under "causal" than under "default"; averaging is one way to get
  # one, since the averaged model inherits the variable order of the *first*
  # model averaged while its edges are the union of all of them.
  data <- test_data()
  data$C <- -data$C
  f1 <- est_DAG(DAG(A ~ B, C ~ C, D ~ D), data, test_tree(), model = "BM")
  f2 <- est_DAG(DAG(C ~ A, B ~ B, D ~ D), data, test_tree(), model = "BM")
  avg <- average_DAGs(list(f1, f2))

  # The set up: acyclic, containing a negative path, and not already sorted.
  expect_true(ggm::isAcyclic((avg$coef != 0) * 1))
  expect_true(any(avg$coef < 0))
  edges <- as.data.frame(avg)
  expect_setequal(paste(edges$from, edges$to), c("B A", "A C"))

  # Ordering by cause puts B before A before C, whatever the matrix order is.
  causal <- as.data.frame(avg, order_by = "causal")
  expect_equal(paste(causal$from, causal$to), c("B A", "A C"))
  # And that is genuinely different from the matrix order, which is what the old
  # code returned here.
  expect_false(identical(rownames(avg$coef), c("B", "A", "C", "D")))
})

test_that("causal ordering gives an informative error for a cyclical model", {
  # Averaging can leave both directions of a path in place when the direction is
  # not resolved, and the result is then not a DAG at all. Previously this
  # surfaced as ggm's bare "The graph is not acyclic!".
  f1 <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  f2 <- est_DAG(DAG(A ~ B, C ~ B), test_data(), test_tree(), model = "BM")
  avg <- average_DAGs(list(f1, f2))
  expect_false(ggm::isAcyclic((avg$coef != 0) * 1))

  expect_error(as.data.frame(avg, order_by = "causal"), "only possible when the model is acyclic")
  # The variables caught in the cycle are named, since that is what the user has
  # to act on. Here A and B point at each other, while C is not involved.
  expect_error(as.data.frame(avg, order_by = "causal"), "paths between A, B form a cycle")
  # And the alternatives are suggested.
  expect_error(as.data.frame(avg, order_by = "causal"), 'order_by = "default"')

  # The other orderings still work on the same model.
  expect_no_error(as.data.frame(avg))
  expect_no_error(as.data.frame(avg, order_by = "strength"))
  # As does coef_plot, which delegates to this method.
  expect_error(coef_plot(avg, error_bar = "se", order_by = "causal"),
               "only possible when the model is acyclic")
})

# coef() and confint() are thin wrappers on as.data.frame(), so what matters is
# that they name the paths identically and stay in the same order, since users
# are invited to cbind() them.

test_that("coef returns the path coefficients as a named vector", {
  fit <- est_DAG(DAG(C ~ A + B, B ~ A), test_data(), test_tree(), model = "BM")
  cf <- coef(fit)

  expect_type(cf, "double")
  expect_null(dim(cf))
  expect_setequal(names(cf), c("A -> B", "A -> C", "B -> C"))
  for (nm in names(cf)) {
    ends <- strsplit(nm, " -> ", fixed = TRUE)[[1]]
    expect_equal(cf[[nm]], fit$coef[ends[1], ends[2]])
  }
  # Absent paths are omitted, as in as.data.frame().
  expect_length(cf, sum(fit$coef != 0))
})

test_that("confint returns the bootstrap bounds, aligned with coef", {
  fit <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                  model = "BM", boot = 20))
  ci <- confint(fit)

  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(rownames(ci), names(coef(fit)))
  expect_true(all(ci[, 1] < coef(fit)))
  expect_true(all(ci[, 2] > coef(fit)))
  # The pairing that the documentation suggests must line up.
  combined <- cbind(coef = coef(fit), confint(fit))
  expect_equal(colnames(combined), c("coef", "2.5 %", "97.5 %"))
})

test_that("confint selects paths with parm", {
  fit <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                  model = "BM", boot = 20))
  one <- confint(fit, parm = "A -> B")
  expect_equal(dim(one), c(1L, 2L))
  expect_equal(rownames(one), "A -> B")
  # A single path still comes back as a matrix, not a bare vector.
  expect_true(is.matrix(confint(fit, parm = 1)))
})

test_that("confint explains that bootstrapping is needed", {
  fit <- est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(), model = "BM")
  expect_error(confint(fit), "fitted without bootstrapping")
  expect_error(confint(fit), "`boot`")
})

test_that("confint refuses a level it cannot deliver", {
  fit <- suppressWarnings(est_DAG(DAG(B ~ A, C ~ B), test_data(), test_tree(),
                                  model = "BM", boot = 20))
  expect_error(confint(fit, level = 0.9), "only the 0.95 bounds are stored")
  expect_error(confint(fit, level = 0.9), "level = 0.9")
  # The default is accepted, whether given or not.
  expect_no_error(confint(fit, level = 0.95))
})
