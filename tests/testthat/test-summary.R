# summary.phylopath() turns the d-sep tables into the model comparison table:
# k, q, C, p, CICc, delta CICc, likelihood and weight.

test_that("summary reports the documented columns and class", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  s <- summary(p)

  expect_s3_class(s, "phylopath_summary")
  expect_s3_class(s, "data.frame")
  expect_named(s, c("model", "k", "q", "C", "p", "CICc", "delta_CICc", "l", "w"))
  expect_setequal(s$model, c("a", "b"))
  # The method itself guards its input.
  expect_error(summary.phylopath(structure(list(), class = "notaphylopath")))
})

test_that("k, q and C are computed from the right quantities", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  s <- summary(p)

  for (mod in names(p$d_sep)) {
    row <- s[s$model == mod, ]
    ds <- p$d_sep[[mod]]
    # k is the number of independence claims.
    expect_equal(row$k, nrow(ds))
    # q is the number of variables plus the number of edges.
    m <- p$model_set[[mod]]
    expect_equal(row$q, nrow(m) + sum(m))
    # C is Fisher's C over the claim p-values.
    expect_equal(row$C, -2 * sum(log(ds$p)))
    expect_equal(row$p, C_p(row$C, row$k))
    expect_equal(row$CICc, CICc(row$C, row$q, nrow(p$data)))
  }
})

test_that("isolates count as variables but not as edges", {
  # define_model_set() encodes an isolate as `D ~ D`, which must not become a
  # self loop in the adjacency matrix, or q would be inflated.
  ms <- define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D))
  expect_true(all(diag(ms$a) == 0))
  expect_equal(sum(ms$a), 2)

  p <- phylo_path(ms, test_data(), test_tree(), model = "BM")
  # 4 variables + 2 edges
  expect_equal(summary(p)$q, 6)
})

test_that("models are ranked by CICc and deltas are relative to the best", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   c = c(C ~ A, B ~ C), .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  s <- summary(p)

  expect_false(is.unsorted(s$CICc))
  expect_equal(s$delta_CICc[1], 0)
  expect_equal(s$delta_CICc, s$CICc - min(s$CICc))
  expect_true(all(s$delta_CICc >= 0))
})

test_that("likelihoods and weights are consistent and normalised", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A),
                                   c = c(C ~ A, B ~ C), .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  s <- summary(p)

  expect_equal(s$l, exp(-0.5 * s$delta_CICc))
  expect_equal(s$w, s$l / sum(s$l))
  expect_equal(sum(s$w), 1)
  expect_equal(which.max(s$w), 1L)
  expect_equal(s$l[1], 1)
})

test_that("summary is unaffected by the order models are supplied in", {
  data <- test_data()
  tree <- test_tree()
  s1 <- summary(phylo_path(
    define_model_set(a = c(B ~ A, C ~ B), b = c(B ~ A, C ~ A), .common = c(D ~ D)),
    data, tree, model = "BM"))
  s2 <- summary(phylo_path(
    define_model_set(b = c(B ~ A, C ~ A), a = c(B ~ A, C ~ B), .common = c(D ~ D)),
    data, tree, model = "BM"))

  expect_equal(s1[order(s1$model), -1], s2[order(s2$model), -1], ignore_attr = TRUE)
})

test_that("print.phylopath_summary returns its input invisibly", {
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
                  test_data(), test_tree(), model = "BM")
  s <- summary(p)
  expect_output(print(s), "CICc")
  out <- capture.output(res <- withVisible(print(s)))
  expect_false(res$visible)
  expect_equal(res$value, s)
})

