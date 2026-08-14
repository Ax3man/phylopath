# define_model_set() assembles a comparable set of causal models: every model
# must span the same variables, with paths absent in a given model simply
# missing and unused variables present as isolates.

test_that("define_model_set returns a named list of DAGs", {
  ms <- define_model_set(one = c(A ~ B), two = c(B ~ A), .common = c(C ~ B))
  expect_type(ms, "list")
  expect_named(ms, c("one", "two"))
  for (m in ms) expect_s3_class(m, "DAG")
})

test_that("all models in a set share the same variables", {
  # Every model must span the union of variables, so that the models are
  # comparable and the d-sep basis sets are over the same variable set.
  expect_equal(unique(vapply(define_model_set(c(A ~ B, B ~ C), c(A ~ C, B ~ C)),
                             nrow, numeric(1))), 3)
  # A variable that appears in only one model still appears in both matrices.
  ms <- define_model_set(one = c(A ~ B, B ~ C), two = c(A ~ C, C ~ D))
  expect_equal(unique(vapply(ms, nrow, numeric(1))), 4)
  for (m in ms) expect_setequal(rownames(m), c("A", "B", "C", "D"))
  # The variable only used by model two is an isolate in model one.
  expect_equal(sum(ms$one["D", ]) + sum(ms$one[, "D"]), 0)
})

test_that(".common paths are added to every model", {
  ms <- define_model_set(one = c(A ~ B), two = c(A ~ C), .common = c(C ~ D))
  for (m in ms) expect_equal(m["D", "C"], 1)
  expect_equal(unique(vapply(ms, nrow, numeric(1))), 4)
})

test_that("isolates can be introduced from either argument", {
  expect_equal(unique(vapply(define_model_set(A ~ B, c(A ~ B, C ~ C)),
                             nrow, numeric(1))), 3)
  expect_equal(unique(vapply(define_model_set(A ~ B, A ~ B, .common = C ~ C),
                             nrow, numeric(1))), 3)
})

test_that("the models a set describes really differ", {
  ms <- define_model_set(one = c(B ~ A), two = c(A ~ B))
  expect_equal(ms$one["A", "B"], 1)
  expect_equal(ms$two["B", "A"], 1)
  expect_false(identical(unclass(ms$one), unclass(ms$two)))
})
