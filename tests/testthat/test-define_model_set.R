context("test-define_model_set")

test_that("classes are correct (list and DAG)", {
  ml <- define_model_set(x = c(A ~ B), .common = c(C ~ D))
  expect_equal(length(ml), 1)
  expect_equal(class(ml), 'list')
  expect_true(inherits(ml[[1]], 'DAG'))
})

test_that("Dimensions of all model matrices is the same.", {
  #Simple two models
  expect_equal(unique(sapply(
    define_model_set(c(A~B, B~C), c(A~C, B~C)), nrow)), 3)
  #Two models with non-overlapping nodes
  expect_equal(unique(sapply(
    define_model_set(c(A~B, B~C), c(A~C, C~D)), nrow)), 4)
  #Use of the .common
  expect_equal(unique(sapply(
    define_model_set(A~B, A~C, .common = C~D), nrow)), 4)
  #Use of isolate
  expect_equal(unique(sapply(
    define_model_set(A~B, c(A~B, C~C)), nrow)), 3)
  #Use of isolate in .common
  expect_equal(unique(sapply(
    define_model_set(A~B, A~B, .common = C~C), nrow)), 3)
})
