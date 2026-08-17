# End-to-end behaviour of phylo_path(). The numeric checks refit the d-sep
# regressions directly with phylolm/phyloglm, so they verify the package's
# plumbing (basis set, deduplication, p-value extraction) rather than just
# freezing whatever it currently returns.

test_that("phylo_path returns the documented structure", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")

  expect_s3_class(p, "phylopath")
  expect_named(p, c("d_sep", "model_set", "data", "tree", "model", "method",
                    "dots", "warnings"))
  expect_named(p$d_sep, c("a", "b"))
  expect_named(p$d_sep$a, c("d_sep", "p", "phylo_par", "model"))
  expect_equal(p$model, "BM")
  expect_equal(p$method, "logistic_MPLE")
  # p-values are p-values.
  for (ds in p$d_sep) {
    expect_true(all(ds$p >= 0 & ds$p <= 1))
  }
})

test_that("every reported d-sep p-value matches a directly fitted phylolm", {
  data <- test_data()
  tree <- test_tree()
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  data, tree, model = "lambda")

  for (mod in names(p$d_sep)) {
    ds <- p$d_sep[[mod]]
    for (i in seq_len(nrow(ds))) {
      fit <- phylolm::phylolm(stats::as.formula(ds$d_sep[i]), data = data,
                              phy = tree, model = "lambda")
      tbl <- stats::coef(phylolm::summary.phylolm(fit))
      expect_equal(ds$p[i], tbl[nrow(tbl), "p.value"])
      expect_equal(ds$phylo_par[i], fit$optpar)
    }
  }
})

test_that("d-sep models shared between candidate models are matched correctly", {
  # phylo_path fits each unique formula once and then maps the results back to
  # the models that need them. A mistake in that mapping would attach the wrong
  # regression, and therefore the wrong p-value, to a model.
  p <- phylo_path(
    define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A), c = c(C ~ A, B ~ C),
                     .common = c(D ~ D)),
    test_data(), test_tree(), model = "BM"
  )
  # There is genuine sharing to get wrong.
  all_stmts <- unlist(lapply(p$d_sep, `[[`, "d_sep"))
  expect_gt(length(all_stmts), length(unique(all_stmts)))

  for (ds in p$d_sep) {
    for (i in seq_len(nrow(ds))) {
      stored <- ds$model[[i]]
      # The stored fit's own formula must be the statement it is filed under.
      expect_equal(
        paste(deparse(stats::formula(stored)), collapse = " "),
        paste(deparse(stats::as.formula(ds$d_sep[i])), collapse = " ")
      )
      expect_equal(ds$p[i], get_p(stored))
    }
  }
})

test_that("the true generating model is preferred over a distinguishable rival", {
  # test_data() was simulated as A -> B -> C with D independent, so the method
  # must rank that model above one that reroutes the chain to A -> C.
  p <- phylo_path(
    define_model_set(true = c(B ~ A, C ~ B), reroute = c(B ~ A, C ~ A),
                     .common = c(D ~ D)),
    test_data(), test_tree(), model = "BM"
  )
  s <- summary(p)
  expect_equal(s$model[1], "true")
  expect_gt(s$w[s$model == "true"], 0.99)
  # The true model is accepted and the rerouted one rejected.
  expect_gt(s$p[s$model == "true"], 0.05)
  expect_lt(s$p[s$model == "reroute"], 0.05)
})

test_that("Markov equivalent models are not distinguished, as expected", {
  # A -> B -> C and C -> B -> A imply exactly the same single independence
  # claim, A _||_ C | B, so no d-sep test can tell them apart. Both must be
  # accepted, with essentially the same C. This is a property of the method, not
  # a shortcoming of the implementation, and it is worth pinning down so that a
  # future change which appears to "improve" discrimination here is recognised
  # as a bug.
  p <- phylo_path(
    define_model_set(forward = c(B ~ A, C ~ B), backward = c(B ~ C, A ~ B),
                     .common = c(D ~ D)),
    test_data(), test_tree(), model = "BM"
  )
  s <- summary(p)
  expect_gt(min(s$p), 0.05)
  expect_equal(s$q[1], s$q[2])
  # The same number of independence claims in both.
  expect_equal(s$k[1], s$k[2])
})

test_that("binary variables are fitted with phyloglm", {
  # cichlids is entirely binary, so every d-sep regression goes through
  # phyloglm rather than phylolm.
  # cichlids has incomplete rows, so dropping them and pruning the tree is
  # expected here.
  p <- suppressMessages(suppressWarnings(phylo_path(
    define_model_set(a = c(C ~ M, P ~ C), b = c(C ~ M, P ~ M), .common = c(G ~ G, D ~ D)),
    cichlids, cichlids_tree
  )))
  expect_true(all(vapply(p$d_sep$a$model, function(m) inherits(m, "phyloglm"), logical(1))))

  # Refit one claim by hand, using the data and tree phylo_path actually used
  # after dropping incomplete rows and pruning.
  ds <- p$d_sep$a
  d2 <- p$data
  f <- stats::as.formula(ds$d_sep[1])
  d2[[all.vars(f)[1]]] <- as.numeric(d2[[all.vars(f)[1]]]) - 1
  fit <- suppressWarnings(phylolm::phyloglm(f, data = d2, phy = p$tree,
                                            method = "logistic_MPLE"))
  tbl <- stats::coef(phylolm::summary.phyloglm(fit))
  expect_equal(ds$p[1], tbl[nrow(tbl), "p.value"])
})

test_that("continuous and binary variables can be mixed in one analysis", {
  data <- test_data()
  data$bin <- factor(rep(c("no", "yes"), 10))
  p <- suppressWarnings(phylo_path(
    define_model_set(a = c(B ~ A, bin ~ B), .common = c(C ~ C, D ~ D)),
    data, test_tree(), model = "BM"
  ))
  classes <- vapply(p$d_sep$a$model, function(m) class(m)[1], character(1))
  expect_true(any(classes == "phylolm"))
  expect_true(any(classes == "phyloglm"))
})

test_that("the stored data and tree are the cleaned ones", {
  # Rows with NAs are dropped and the tree pruned, and what is stored has to be
  # the reduced data, since everything downstream (n in CICc, est_DAG) uses it.
  data <- test_data()
  data$A[1:2] <- NA
  p <- suppressMessages(phylo_path(
    define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
    data, test_tree(), model = "BM"
  ))
  expect_equal(nrow(p$data), 38)
  expect_length(p$tree$tip.label, 38)
  expect_false("t1" %in% p$tree$tip.label)
  # Unused columns are gone too.
  expect_setequal(names(p$data), c("A", "B", "C", "D"))
})

test_that("na.rm = FALSE refuses incomplete data", {
  data <- test_data()
  data$A[1] <- NA
  expect_error(
    phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
               data, test_tree(), model = "BM", na.rm = FALSE),
    "NA values were found"
  )
})

test_that("a supplied `order` is used and validated", {
  data <- test_data()
  tree <- test_tree()
  ms <- define_model_set(a = c(C ~ B), .common = c(A ~ A, D ~ D))

  # A and D are isolates, so the order alone decides which is the response.
  p1 <- phylo_path(ms, data, tree, model = "BM", order = c("A", "D", "B", "C"))
  p2 <- phylo_path(ms, data, tree, model = "BM", order = c("D", "A", "B", "C"))
  stmt1 <- grep("^D ~ |^A ~ ", p1$d_sep$a$d_sep, value = TRUE)
  stmt2 <- grep("^D ~ |^A ~ ", p2$d_sep$a$d_sep, value = TRUE)
  expect_true(any(grepl("^D ~ .*A$", stmt1)))
  expect_true(any(grepl("^A ~ .*D$", stmt2)))

  expect_error(
    phylo_path(ms, data, tree, model = "BM", order = c("A", "A", "B", "C")),
    "must contain each variable exactly once"
  )
  # An incomplete order, and one naming variables that are not in the models,
  # both have to be reported by name rather than failing somewhere downstream.
  expect_error(
    phylo_path(ms, data, tree, model = "BM", order = c("A", "B", "C")),
    "Missing from `order`: D"
  )
  expect_error(
    phylo_path(ms, data, tree, model = "BM", order = c("X1", "X2", "X3", "X4")),
    "Not a variable in the models: X1, X2, X3, X4"
  )
  expect_error(
    phylo_path(ms, data, tree, model = "BM", order = c("A", "B", "C", "typo")),
    "Missing from `order`: D"
  )
})

test_that("superseded and deprecated arguments warn", {
  ms <- define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D))
  expect_warning(
    phylo_path(ms, test_data(), test_tree(), model = "BM", boot = 10),
    "`boot` has been deprecated here"
  )
  expect_warning(
    phylo_path(ms, test_data(), test_tree(), model = "BM", parallel = 2),
    "superseded by the use of `future`"
  )
  # Passing boot must not leak it into the fitted models.
  p <- suppressWarnings(phylo_path(ms, test_data(), test_tree(), model = "BM", boot = 10))
  expect_length(p$dots, 0)
})

# Collect every warning a call emits, rather than just the first, so that the
# *number* of warnings can be asserted.
all_warnings <- function(expr) {
  ws <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  ws
}

test_that("unrecognised arguments in ... warn once, and name themselves", {
  ms <- define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D))
  ws <- all_warnings(
    p <- phylo_path(ms, test_data(), test_tree(), model = "BM", not_an_argument = 1)
  )
  # Once per user call, not once per regression: this model set has four.
  expect_length(ws, 1)
  expect_match(ws, "not recognized", fixed = FALSE)
  # The message says which argument was wrong, which is the actionable part.
  expect_match(ws, "not_an_argument")
  # The offending argument is discarded, so it cannot reach phylolm later.
  expect_length(p$dots, 0)
})

test_that("unrecognised arguments warn once from the fitting functions too", {
  ms <- define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A), .common = c(D ~ D))
  p <- phylo_path(ms, test_data(), test_tree(), model = "BM")
  d <- DAG(B ~ A, C ~ B)

  expect_length(all_warnings(est_DAG(d, test_data(), test_tree(), model = "BM",
                                     nonsense = 1)), 1)
  expect_length(all_warnings(best(p, nonsense = 1)), 1)
  expect_length(all_warnings(choice(p, "a", nonsense = 1)), 1)
  # average() fits one DAG per model within the cut off, so this is the case that
  # used to warn once per model per regression. The cut off is set explicitly so
  # that more than one model is certain to be averaged over.
  expect_length(all_warnings(average(p, cut_off = Inf, nonsense = 1)), 1)
})

test_that("arguments recognised by only one of phylolm and phyloglm are kept", {
  # `btol` is a phyloglm argument and `measurement_error` a phylolm one. A data
  # set with both binary and continuous variables needs both to survive the
  # check, since narrowing to the relevant one happens per fit.
  ms <- define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D))
  ws <- all_warnings(
    p <- phylo_path(ms, test_data(), test_tree(), model = "BM",
                    btol = 20, measurement_error = FALSE)
  )
  expect_length(ws, 0)
  expect_named(p$dots, c("btol", "measurement_error"))
})

test_that("errors from the underlying model are reported with the formula", {
  expect_error(
    phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
               test_data(), test_tree(), model = "not_a_model"),
    "Fitting the following model"
  )
})

test_that("arguments are passed through to phylolm", {
  ms <- define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D))
  bm <- phylo_path(ms, test_data(), test_tree(), model = "BM")
  lambda <- phylo_path(ms, test_data(), test_tree(), model = "lambda")
  # Brownian motion has no free phylogenetic parameter; Pagel's lambda does.
  expect_true(all(is.na(bm$d_sep$a$phylo_par)))
  expect_false(any(is.na(lambda$d_sep$a$phylo_par)))
  expect_false(isTRUE(all.equal(bm$d_sep$a$p, lambda$d_sep$a$p)))
})

test_that("print.phylopath reports the variables and models", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  out <- capture.output(print(p))
  expect_match(paste(out, collapse = " "), "Continuous:")
  for (v in c("A", "B", "C", "D")) expect_match(paste(out, collapse = " "), v)
  expect_match(paste(out, collapse = " "), "these models: a b")
  expect_match(paste(out, collapse = " "), "phylogenetic regressions")
})

test_that("a single model warning is still announced", {
  # The notice has to fire for one warning as well as for several, or the user
  # never learns that show_warnings() has something to show.
  data <- test_data()
  data$b1 <- factor(rep(c("no", "yes"), 20))
  data$b2 <- factor(rep(c("no", "yes"), each = 20))
  ms <- define_model_set(a = c(B ~ A), .common = c(C ~ C, D ~ D, b1 ~ b1, b2 ~ b2))

  p <- suppressWarnings(phylo_path(ms, data, test_tree(), model = "BM"))
  expect_length(p$warnings, 1)
  expect_warning(phylo_path(ms, data, test_tree(), model = "BM"),
                 "Use `show_warnings\\(\\)` to view them")
})

