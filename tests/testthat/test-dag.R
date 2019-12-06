context("test-dag")
library(phylopath)

test_that("DAGs get the correct classes", {
  expect_true(inherits(DAG(A~B, B~C), 'matrix'))
  expect_true(inherits(DAG(A~B, B~C), 'DAG'))
})

test_that('DAGs make correct nr of rows and columns', {
  expect_equal(nrow(DAG(A~B, B~C, C~D)), 4)
  expect_equal(nrow(DAG(A~B, B~C, C~D)), ncol(DAG(A~B, B~C, C~D)))
})

test_that('Correct ordering of DAGs', {
  expect_equal(
    unclass(DAG(A~B, B~C, C~D)),
    structure(
      c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
      .Dim = c(4L, 4L),
      .Dimnames = list(c("D", "C", "B", "A"), c("D", "C", "B", "A"))
    )
  )
  expect_equal(
    unclass(DAG(A~B, C~D, B~C, order = FALSE)),
    structure(
      c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
      .Dim = c(4L, 4L), .Dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
    )
  )
})
