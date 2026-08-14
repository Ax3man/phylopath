# find_consensus_order() decides the causal order used to orient d-sep
# regressions for variable pairs that no model links causally. Because PGLS
# p-values are not symmetric in the two variables, this choice affects results,
# so the order has to be a genuine topological order whenever one exists.

# Does every edge of `m` point forwards in `order`?
edges_go_forward <- function(m, order) {
  idx <- stats::setNames(seq_along(order), order)
  for (from in rownames(m)) {
    for (to in colnames(m)) {
      if (m[from, to] == 1 && from != to && idx[[from]] >= idx[[to]]) return(FALSE)
    }
  }
  TRUE
}

test_that("the union of the models is used when it is acyclic", {
  # A -> B -> C and A -> C combine into an acyclic graph, whose unique
  # topological order is A, B, C.
  ms <- define_model_set(one = c(B ~ A, C ~ B), two = c(B ~ A, C ~ A))
  expect_equal(find_consensus_order(ms), c("A", "B", "C"))

  ms2 <- define_model_set(one = c(X2 ~ X1, X3 ~ X2), two = c(X3 ~ X2, X4 ~ X3),
                          .common = c(X2 ~ X1))
  ord <- find_consensus_order(ms2)
  expect_setequal(ord, c("X1", "X2", "X3", "X4"))
  for (m in ms2) expect_true(edges_go_forward(m, ord))
})

test_that("the order is a topological sort of every model when the union is acyclic", {
  set.seed(11)
  checked <- 0
  for (i in 1:25) {
    v <- sample(3:5, 1)
    ms <- replicate(sample(2:4, 1), random_DAG(v), simplify = FALSE)
    same <- lapply(ms, function(x) x[rownames(ms[[1]]), colnames(ms[[1]])])
    if (!ggm::isAcyclic(sign(Reduce(`+`, same)))) next
    ord <- find_consensus_order(ms)
    for (m in ms) expect_true(edges_go_forward(m, ord))
    checked <- checked + 1
  }
  expect_gt(checked, 0)
})

test_that("a cyclic union falls back to the majority ordering", {
  # B ~ A and A ~ B cannot both be respected; the result must still be a usable
  # total order over the variables.
  ms <- define_model_set(one = c(B ~ A, C ~ B), two = c(A ~ B, C ~ B))
  ord <- find_consensus_order(ms)
  expect_setequal(ord, c("A", "B", "C"))
  expect_length(ord, 3)

  # With three models favouring A -> B and one favouring B -> A, A wins.
  ms2 <- define_model_set(
    a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A), c = c(B ~ A, A ~ C), d = c(A ~ B, C ~ B)
  )
  ord2 <- find_consensus_order(ms2)
  expect_lt(which(ord2 == "A"), which(ord2 == "B"))
})

test_that("the order is always a permutation of the variables", {
  # This also exercises the "contact the maintainer (code 2)" branch, which
  # should be unreachable.
  set.seed(12)
  for (i in 1:25) {
    v <- sample(3:5, 1)
    ms <- replicate(sample(2:4, 1), random_DAG(v, p = stats::runif(1, 0.2, 0.7)),
                    simplify = FALSE)
    ord <- find_consensus_order(ms)
    expect_setequal(ord, colnames(ms[[1]]))
    expect_length(ord, v)
    expect_null(names(ord))
  }
})

test_that("models whose rows and columns are in different orders are handled", {
  # define_model_set() topologically sorts each model separately, so the
  # matrices in a model set routinely disagree on row order.
  ms <- define_model_set(one = c(B ~ A, C ~ B), two = c(A ~ B, C ~ A))
  expect_false(identical(rownames(ms[[1]]), rownames(ms[[2]])))
  expect_setequal(find_consensus_order(ms), c("A", "B", "C"))
})

test_that("a single model gives its own topological order", {
  ms <- define_model_set(only = c(B ~ A, C ~ B, D ~ C))
  expect_equal(find_consensus_order(ms), c("A", "B", "C", "D"))
})
