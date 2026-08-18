# The information-theoretic machinery of the method, checked against the
# published equations in von Hardenberg & Gonzalez-Voyer (2013),
# doi:10.1111/j.1558-5646.2012.01790.x

test_that("C_stat implements Fisher's C = -2 * sum(ln(p))", {
  expect_equal(C_stat(c(0.5, 0.5)), -2 * (log(0.5) + log(0.5)))
  expect_equal(C_stat(0.05), -2 * log(0.05))
  # A perfectly fitting model (all independencies hold exactly) gives C = 0.
  expect_equal(C_stat(c(1, 1, 1)), 0)
  # C grows as the independence claims are violated more strongly.
  expect_gt(C_stat(c(0.01, 0.01)), C_stat(c(0.4, 0.4)))
  # An empty basis set (saturated model) contributes nothing.
  expect_equal(C_stat(numeric(0)), 0)
})

test_that("C is tested against a chi-square with 2k degrees of freedom", {
  expect_equal(C_p(5, 3), stats::pchisq(5, df = 6, lower.tail = FALSE))
  expect_equal(C_p(0, 2), 1)
  # A larger C at fixed k is less probable.
  expect_lt(C_p(20, 3), C_p(5, 3))
  # More claims at fixed C is more probable, because df grows.
  expect_gt(C_p(10, 6), C_p(10, 2))
})

test_that("CICc implements C + 2q(n / (n - q - 1))", {
  expect_equal(CICc(5, 4, 50), 5 + 2 * 4 * (50 / (50 - 4 - 1)))
  # With no free parameters CICc reduces to C.
  expect_equal(CICc(7, 0, 50), 7)
  # The small sample penalty grows as n shrinks towards q.
  expect_gt(CICc(5, 6, 20) - 5, CICc(5, 6, 200) - 5)
})

test_that("likelihoods and weights follow the standard IT formulae", {
  expect_equal(l(0), 1)
  expect_equal(l(4), exp(-2))
  expect_equal(w(c(1, 1)), c(0.5, 0.5))
  expect_equal(sum(w(l(c(0, 1, 5, 10)))), 1)
  # The best model (delta = 0) always carries the greatest weight.
  ws <- w(l(c(0, 2, 4)))
  expect_equal(which.max(ws), 1L)
})

test_that("published Table 2 of von Hardenberg & Gonzalez-Voyer is reproduced", {
  # Table 2a, the 68 species avian broodmate competition analysis. C is reported
  # to two decimals in the paper, which bounds how exactly we can reproduce the
  # derived quantities.
  #
  # Models O and A are omitted: the paper lists q = 9 for both, but its CICc
  # column for those two rows is the value you get with q = 10. Every other row
  # is internally consistent, and k + q = 15 throughout (as it must be, since
  # k = choose(5, 2) - edges and q = 5 + edges), which confirms q = 9 is the
  # correct value for them and the printed CICc is the typo.
  paper <- data.frame(
    model = c("K", "I", "B", "L", "M", "C", "D", "F", "N", "G", "E", "H"),
    C     = c(5.28, 8.83, 9.37, 10.01, 11.18, 14.11, 12.26, 13.02, 16.64, 17.11, 19.04, 26.65),
    k     = c(4, 5, 5, 5, 5, 6, 5, 5, 6, 6, 6, 6),
    q     = c(11, 10, 10, 10, 10, 9, 10, 10, 9, 9, 9, 9),
    p     = c(0.727, 0.548, 0.497, 0.439, 0.343, 0.294, 0.268, 0.222, 0.163, 0.145, 0.087, 0.008),
    CICc  = c(31.994, 32.693, 33.230, 33.880, 35.043, 35.212, 36.122, 36.880, 37.748, 38.212,
              40.145, 47.749)
  )
  n <- 68

  expect_equal(C_p(paper$C, paper$k), paper$p, tolerance = 0.002)
  expect_equal(CICc(paper$C, paper$q, n), paper$CICc, tolerance = 0.02)

  # k and q are not independent: for a model set on v variables,
  # k = choose(v, 2) - edges while q = v + edges, so k + q = choose(v, 2) + v.
  expect_true(all(paper$k + paper$q == choose(5, 2) + 5))
})

test_that("par_avg reproduces the weighted average and unconditional SE", {
  # Reference values computed by hand from the MuMIn::par.avg formulation:
  # the averaged SE combines within-model variance and between-model spread.
  x <- c(1, 2)
  se <- c(0.1, 0.3)
  wt <- c(0.75, 0.25)
  res <- par_avg(x, se, wt)

  expect_named(res, c("Coefficient", "SE", "Lower CI", "Upper CI"))
  expect_equal(unname(res["Coefficient"]), 1.25)
  expect_equal(
    unname(res["SE"]),
    sqrt(0.75 * (0.1^2 + (1 - 1.25)^2) + 0.25 * (0.3^2 + (2 - 1.25)^2))
  )
  # Intervals are symmetric normal intervals around the averaged coefficient.
  expect_equal(unname(res["Lower CI"]), 1.25 - stats::qnorm(0.975) * unname(res["SE"]))
  expect_equal(unname(res["Upper CI"]), 1.25 + stats::qnorm(0.975) * unname(res["SE"]))
})

test_that("par_avg ignores NA coefficients and renormalises the weights", {
  # This is what conditional averaging relies on: a path absent from one model
  # is NA there, and must be averaged over the remaining models only.
  both <- par_avg(c(2, 2), c(0.2, 0.2), c(0.9, 0.1))
  one  <- par_avg(c(2, NA), c(0.2, NA), c(0.9, 0.1))
  expect_equal(unname(one["Coefficient"]), 2)
  expect_equal(unname(one["SE"]), 0.2)
  expect_equal(both[["Coefficient"]], one[["Coefficient"]])
})

test_that("par_avg validates its inputs", {
  expect_error(par_avg("a", 1, 1), "must be numeric vectors")
  expect_error(par_avg(c(1, 2), 1, c(1, 1)), "must all be the same length")
  expect_error(par_avg(c(1, 2), c(1, 1), 1), "must all be the same length")
})

test_that("C_p keeps its precision far into the upper tail", {
  # Computing the tail as 1 - pchisq() underflows to exactly 0 once the tail is
  # below about 1e-16, which reports an impossible p-value of 0 for a model that
  # is merely very strongly rejected.
  expect_gt(C_p(200, 5), 0)
  expect_equal(C_p(200, 5), stats::pchisq(200, 10, lower.tail = FALSE))
  expect_equal(1 - stats::pchisq(200, 10), 0)

  # Still monotonic, and unchanged where the naive form was accurate anyway.
  expect_lt(C_p(300, 5), C_p(200, 5))
  expect_equal(C_p(5, 3), 1 - stats::pchisq(5, 6))
  expect_equal(C_p(0, 2), 1)
})

test_that("the p-value floor keeps Fisher's C finite", {
  # get_p() floors each claim at the smallest representable number so that
  # C = -2 * sum(log(p)) cannot become infinite. An infinite C would make CICc
  # infinite, and then every model weight NaN.
  expect_false(is.finite(C_stat(0)))
  expect_true(is.finite(C_stat(.Machine$double.xmin)))
  expect_true(is.finite(C_stat(rep(.Machine$double.xmin, 10))))
  expect_equal(C_stat(.Machine$double.xmin), -2 * log(.Machine$double.xmin))
})
