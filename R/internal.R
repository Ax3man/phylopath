set_to_formula <- function(x) {
  stats::formula(paste(x[1], paste(x[-1], collapse = '+'), sep = '~'))
}

find_formulas <- function(d) {
  s <- ggm::basiSet(d)
  lapply(s, set_to_formula)
}

C_stat <- function(ps) -2 * sum(log(ps))

C_p <- function(C, k) 1 - stats::pchisq(C, 2 * k)

get_p <- function(m) summary(m)$tTable[2, 'p-value']

CICc <- function(C, q, n) C + 2 * q * (n / (n - 1 - q))

l <- function(dCICc) exp(-0.5 * dCICc)

w <- function(l) l / sum(l)

gls2 <- function(..., cor_fun, tree) {
  for (par in seq(1, 0, -0.05)) {
    m <- tryCatch(
      nlme::gls(..., correlation = cor_fun(par, tree)),
      error = function(e) NA, warning = function(w) NA)
    if ('gls' %in% class(m)) break
  }
  return(m)
}