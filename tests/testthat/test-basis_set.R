# The d-separation basis set: the heart of the method. A wrong basis set, or a
# regression equation with the wrong response or the tested variable in the
# wrong position, would silently produce wrong p-values and so wrong C, CICc
# and model ranking.

# Helper: the pair being tested by a formula is (response, last term). The last
# term matters because get_p() reads the p-value off the *last* coefficient.
formula_pair <- function(f) {
  v <- all.vars(f)
  c(v[1], v[length(v)])
}
formula_cond <- function(f) {
  v <- all.vars(f)
  v[-c(1, length(v))]
}

test_that("the basis set reproduces Figure 1A of von Hardenberg & Gonzalez-Voyer", {
  # The paper spells out the basis set for the chain X1 -> X2 -> X3 -> X4 -> X5:
  #   (X1,X3){X2}  (X1,X4){X3}  (X1,X5){X4}
  #   (X2,X4){X1,X3}  (X2,X5){X1,X4}  (X3,X5){X2,X4}
  # A chain has a unique topological order, so the regressions are fully
  # determined: the downstream member of each pair is the response, and the
  # tested member comes last.
  d <- DAG(X2 ~ X1, X3 ~ X2, X4 ~ X3, X5 ~ X4)
  expect_equal(rownames(d), c("X1", "X2", "X3", "X4", "X5"))

  formulas <- find_formulas(d, rownames(d))
  expect_equal(
    vapply(formulas, function(f) paste(deparse(f), collapse = " "), character(1)),
    c("X3 ~ X2 + X1", "X4 ~ X3 + X1", "X5 ~ X4 + X1",
      "X4 ~ X1 + X3 + X2", "X5 ~ X1 + X4 + X2", "X5 ~ X2 + X4 + X3")
  )
})

test_that("the basis set has choose(v, 2) - edges claims", {
  # Every pair of variables is either directly linked, or generates exactly one
  # independence claim. This is the identity that makes k and q complementary
  # (see test-statistics.R).
  expect_length(find_formulas(DAG(X2 ~ X1, X3 ~ X2, X4 ~ X3, X5 ~ X4), LETTERS), choose(5, 2) - 4)
  expect_length(find_formulas(DAG(B ~ A, C ~ B), c("A", "B", "C")), choose(3, 2) - 2)

  set.seed(1)
  for (i in 1:15) {
    v <- sample(4:6, 1)
    d <- random_DAG(v)
    n_edges <- sum(d)
    # A saturated DAG has no claims to test and is rejected elsewhere.
    if (n_edges == choose(v, 2)) next
    expect_length(find_formulas(d, rownames(d)), choose(v, 2) - n_edges)
  }
})

test_that("each claim conditions on the parents of both variables", {
  # "one lists the parent variables of either nonadjacent variables"
  set.seed(2)
  for (i in 1:15) {
    d <- random_DAG(sample(4:6, 1))
    if (sum(d) == choose(ncol(d), 2)) next
    for (f in find_formulas(d, rownames(d))) {
      pair <- formula_pair(f)
      parents <- rownames(d)[d[, pair[1]] == 1 | d[, pair[2]] == 1]
      expect_setequal(formula_cond(f), setdiff(parents, pair))
    }
  }
})

test_that("the response is the downstream variable whenever a path exists", {
  # If X causes Y (directly or indirectly) the claim must be tested by
  # regressing Y on X, never the other way round.
  set.seed(3)
  for (i in 1:15) {
    d <- random_DAG(sample(4:6, 1))
    if (sum(d) == choose(ncol(d), 2)) next
    for (f in find_formulas(d, rownames(d))) {
      pair <- formula_pair(f)
      resp <- pair[1]
      tested <- pair[2]
      # There must be no directed path from the response to the tested variable.
      path_from_response <- ggm::findPath(
        d, which(rownames(d) == resp), which(rownames(d) == tested)
      )
      expect_null(path_from_response)
    }
  }
})

test_that("the tested variable is last, which is where get_p() looks", {
  # get_p() reads the p-value from the final row of the coefficient table, so
  # the tested variable must be the final term of the formula.
  d <- DAG(X2 ~ X1, X3 ~ X2, X4 ~ X3, X5 ~ X4)
  for (f in find_formulas(d, rownames(d))) {
    terms <- all.vars(f)
    expect_false(terms[length(terms)] %in% terms[-length(terms)])
  }

  # Confirm the linkage end to end: the reported p-value equals the one for the
  # last coefficient of a directly fitted model.
  tree <- test_tree()
  data <- test_data()
  p <- phylo_path(define_model_set(a = c(B ~ A, C ~ B), .common = c(D ~ D)),
                  data, tree, model = "BM")
  ds <- p$d_sep$a
  for (i in seq_len(nrow(ds))) {
    f <- stats::as.formula(ds$d_sep[i])
    fit <- phylolm::phylolm(f, data = data, phy = tree, model = "BM")
    tbl <- stats::coef(phylolm::summary.phylolm(fit))
    expect_equal(ds$p[i], tbl[nrow(tbl), "p.value"])
    # and that really is the coefficient of the tested variable
    expect_equal(rownames(tbl)[nrow(tbl)], formula_pair(f)[2])
  }
})

test_that("`order` decides the direction for pairs with no causal path", {
  # A and B are both isolates here, so nothing in the graph orients the claim
  # and the supplied causal order has to.
  d <- DAG(C ~ D, A ~ A, B ~ B)
  f_ab <- Filter(
    function(f) setequal(formula_pair(f), c("A", "B")),
    find_formulas(d, c("A", "B", "C", "D"))
  )
  expect_length(f_ab, 1)
  expect_equal(formula_pair(f_ab[[1]]), c("B", "A"))

  f_ba <- Filter(
    function(f) setequal(formula_pair(f), c("A", "B")),
    find_formulas(d, c("B", "A", "C", "D"))
  )
  expect_equal(formula_pair(f_ba[[1]]), c("A", "B"))
})

test_that("set_to_formula puts the response first and the tested variable last", {
  expect_equal(
    paste(deparse(set_to_formula(c("X", "Y", "Z1", "Z2"))), collapse = " "),
    "Y ~ Z1 + Z2 + X"
  )
  expect_equal(paste(deparse(set_to_formula(c("X", "Y"))), collapse = " "), "Y ~ X")
})

test_that("a saturated model cannot be tested", {
  # Every pair is adjacent, so the basis set is empty and there is nothing to
  # test. This must be an informative error, not a silent C = 0.
  expect_error(
    find_formulas(DAG(B ~ A, C ~ A + B), c("A", "B", "C")),
    "fully connected"
  )
  expect_error(
    phylo_path(define_model_set(a = c(B ~ A, C ~ A + B)), test_data(), test_tree()),
    "fully connected"
  )
})

test_that("find_formulas never hits its internal consistency errors", {
  # find_formulas has two "contact the maintainer" branches, guarding against
  # basiSet returning a reversed pair and against cycles. Neither should ever
  # trigger for a valid DAG.
  set.seed(4)
  for (i in 1:40) {
    d <- random_DAG(sample(3:6, 1), p = stats::runif(1, 0.2, 0.7))
    if (sum(d) == choose(ncol(d), 2)) next
    expect_no_error(find_formulas(d, rownames(d)))
  }
})
