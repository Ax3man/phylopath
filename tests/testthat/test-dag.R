# DAG() builds the adjacency matrices everything else depends on. Row/column
# order carries meaning (it is the topological order), so it is tested
# explicitly.

test_that("DAG returns a matrix with the DAG class", {
  d <- DAG(A ~ B, B ~ C)
  expect_true(is.matrix(d))
  expect_s3_class(d, "DAG")
  expect_true(all(d %in% c(0, 1)))
})

test_that("DAG is square with one row per variable", {
  d <- DAG(A ~ B, B ~ C, C ~ D)
  expect_equal(dim(d), c(4, 4))
  expect_equal(rownames(d), colnames(d))
})

test_that("edges point from parent (row) to child (column)", {
  d <- DAG(B ~ A)
  expect_equal(d["A", "B"], 1)
  expect_equal(d["B", "A"], 0)
  # Multiple parents in one formula.
  d2 <- DAG(C ~ A + B)
  expect_equal(d2["A", "C"], 1)
  expect_equal(d2["B", "C"], 1)
  expect_equal(sum(d2), 2)
})

test_that("nodes are in topological order by default", {
  expect_equal(rownames(DAG(A ~ B, B ~ C, C ~ D)), c("D", "C", "B", "A"))
  # With order = FALSE the order of first appearance is kept.
  expect_equal(rownames(DAG(A ~ B, C ~ D, B ~ C, order = FALSE)),
               c("A", "B", "C", "D"))
  # Topological order means every edge points forward.
  expect_true(edges_go_forward(DAG(A ~ B, B ~ C, C ~ D)))
})

test_that("a variable declared as its own parent is an isolate, not a self loop", {
  d <- DAG(A ~ B + C, D ~ D)
  expect_true("D" %in% rownames(d))
  expect_equal(sum(d["D", ]), 0)
  expect_equal(sum(d[, "D"]), 0)
  expect_true(all(diag(d) == 0))
})

test_that("a cyclic specification is rejected", {
  # ggm::DAG refuses to topologically sort a graph with a cycle.
  expect_error(suppressWarnings(DAG(A ~ B, B ~ C, C ~ A)))
  expect_warning(try(DAG(A ~ B, B ~ C, C ~ A), silent = TRUE), "directed cycles")
})
