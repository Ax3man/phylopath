# Input validation. These checks are the package's first line of defence, and
# most of them exist to turn confusing downstream failures in phylolm into
# actionable messages.

test_that("all models must contain the same variables", {
  ms <- list(DAG(B ~ A), DAG(C ~ A))
  d <- tiny_data(A = 1:4, B = 5:8, C = 9:12)
  expect_error(
    check_models_data_tree(ms, d, tiny_tree(), na.rm = FALSE),
    "All causal models need to include the same variables"
  )
  # The message lists the offending variables, to help the user spot the typo.
  expect_error(
    check_models_data_tree(ms, d, tiny_tree(), na.rm = FALSE),
    "A\nB\nC"
  )
})

test_that("a single model is never rejected for inconsistent variables", {
  res <- check_models_data_tree(list(DAG(B ~ A)), tiny_data(A = 1:4, B = 5:8),
                               tiny_tree(), na.rm = FALSE)
  expect_named(res, c("model_set", "data", "tree"))
})

test_that("columns not used by any model are dropped", {
  d <- tiny_data(A = 1:4, B = 5:8, unused = 9:12)
  res <- check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE)
  expect_named(res$data, c("A", "B"))
  expect_equal(rownames(res$data), c("a", "b", "c", "d"))
})

test_that("character columns become factors", {
  d <- tiny_data(A = c("x", "y", "x", "y"), B = c("p", "p", "q", "q"))
  res <- check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE)
  expect_true(is.factor(res$data$A))
  expect_true(is.factor(res$data$B))
  # Level order is preserved, which fixes which level phyloglm treats as 1.
  expect_equal(levels(res$data$A), c("x", "y"))
})

test_that("binary factors are accepted and other factors are rejected", {
  two <- tiny_data(A = 1:4, B = factor(c("x", "y", "x", "y")))
  expect_no_error(check_models_data_tree(list(DAG(B ~ A)), two, tiny_tree(), na.rm = FALSE))

  three <- tiny_data(A = 1:4, B = factor(c("x", "y", "z", "x")))
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), three, tiny_tree(), na.rm = FALSE),
    "Variable 'B' is expected to b.*inary, but has 3 levels"
  )
  one <- tiny_data(A = 1:4, B = factor(rep("x", 4)))
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), one, tiny_tree(), na.rm = FALSE),
    "has 1 levels"
  )
})

test_that("numeric 0/1 columns are flagged as probable user error", {
  # A binary variable passed as numeric would silently be modelled as
  # continuous, so the user is warned to convert it to a factor.
  d <- tiny_data(A = c(0, 1, 0, 1), B = c(2.5, 3.5, 4.5, 5.5))
  expect_warning(
    check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE),
    "Column A appears to have binary data"
  )
  # Genuinely continuous data draws no warning.
  expect_no_warning(
    check_models_data_tree(list(DAG(B ~ A)), tiny_data(A = c(0.5, 1, 2, 3), B = 5:8),
                           tiny_tree(), na.rm = FALSE)
  )
})

test_that("the tree must be a single phylo object", {
  d <- tiny_data(A = 1:4, B = 5:8)
  multi <- structure(list(tiny_tree(), tiny_tree()), class = "multiPhylo")
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), d, multi, na.rm = FALSE),
    "You are passing several trees"
  )
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), d, "not a tree", na.rm = FALSE),
    "tree needs to be of class `phylo`"
  )
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), d, NULL, na.rm = FALSE),
    "tree needs to be of class `phylo`"
  )
})

test_that("na.rm controls whether incomplete rows are dropped or refused", {
  d <- tiny_data(A = c(1, NA, 3, 4), B = 5:8)
  # Two messages: the dropped row, and the tree being pruned to match.
  expect_message(
    expect_message(
      res <- check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = TRUE),
      "1 rows were dropped"
    ),
    "Pruned tree"
  )
  expect_equal(nrow(res$data), 3)
  expect_equal(rownames(res$data), c("a", "c", "d"))
  # The tree is pruned to match the reduced data.
  expect_setequal(res$tree$tip.label, c("a", "c", "d"))

  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE),
    "NA values were found"
  )
})

test_that("NAs only matter in the variables actually used", {
  # The unused column is dropped before the NA check, so it must not trigger it.
  d <- tiny_data(A = 1:4, B = 5:8, unused = c(NA, NA, NA, NA))
  expect_no_error(
    check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE)
  )
})

test_that("data rownames must be matched by tip labels", {
  d <- data.frame(A = 1:4, B = 5:8, row.names = c("w", "x", "y", "z"))
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE),
    "rownames that are exactly matched by name with tips in the tree"
  )
  # A species missing from the tree is an error, not a silent drop.
  d2 <- data.frame(A = 1:5, B = 5:9, row.names = c("a", "b", "c", "d", "e"))
  expect_error(
    check_models_data_tree(list(DAG(B ~ A)), d2, tiny_tree(), na.rm = FALSE),
    "rownames that are exactly matched"
  )
})

test_that("tips absent from the data are pruned with a message", {
  d <- data.frame(A = 1:2, B = 3:4, row.names = c("a", "c"))
  expect_message(
    res <- check_models_data_tree(list(DAG(B ~ A)), d, tiny_tree(), na.rm = FALSE),
    "Pruned tree"
  )
  expect_setequal(res$tree$tip.label, c("a", "c"))
  # A tree that already matches is left alone and stays quiet.
  full <- tiny_data(A = 1:4, B = 5:8)
  expect_no_message(
    check_models_data_tree(list(DAG(B ~ A)), full, tiny_tree(), na.rm = FALSE)
  )
})

test_that("unnamed model sets are lettered", {
  ms <- list(DAG(B ~ A), DAG(B ~ A))
  res <- check_models_data_tree(ms, tiny_data(A = 1:4, B = 5:8), tiny_tree(), na.rm = FALSE)
  expect_equal(names(res$model_set), c("A", "B"))

  # Existing names are left untouched.
  named <- list(first = DAG(B ~ A), second = DAG(B ~ A))
  res2 <- check_models_data_tree(named, tiny_data(A = 1:4, B = 5:8), tiny_tree(), na.rm = FALSE)
  expect_equal(names(res2$model_set), c("first", "second"))
})

test_that("a variable used in the models but absent from the data is named", {
  # Previously this surfaced as "undefined columns selected" from [.data.frame.
  d <- tiny_data(A = 1:4, B = 5:8)
  expect_error(
    check_models_data_tree(list(DAG(B ~ A, C ~ A)), d, tiny_tree(), na.rm = FALSE),
    "are used in the causal models, but are not columns in `data`: C"
  )
  # Several missing variables are all listed (in the DAG's own node order).
  msg <- tryCatch(
    check_models_data_tree(list(DAG(B ~ A, C ~ A, E ~ A)), d, tiny_tree(), na.rm = FALSE),
    error = conditionMessage
  )
  expect_match(msg, "C")
  expect_match(msg, "E")
  expect_match(msg, ", ")
})

