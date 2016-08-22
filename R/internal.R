find_consensus_order <- function(models) {
  # If the fully combined model is acyclic, then we use that.
  full_model <- sign(Reduce('+', models))
  if (ggm::isAcyclic(full_model)) {
    return(rownames(ggm::topSort(full_model)))
  }
  # Otherwise we find the most common orderings and use those.
  vars <- lapply(models, row.names)
  combs <- as.data.frame(t(utils::combn(vars[[1]], 2)), stringsAsFactors = FALSE)
  combs$count <- 0
  for (i in seq_along(vars)) {
    v <- apply(combs, 1, function(x) {
      which(vars[[i]] == x[1]) < which(vars[[i]] == x[2])
    } )
    combs$count <- combs$count + v
  }
  combs <- dplyr::mutate_(combs,
                          V1 = ~ifelse(count > q, V1, V2),
                          V2 = ~ifelse(count > q, V2, V1))
  combs <- dplyr::group_by_(combs, ~V1)
  combs <- dplyr::mutate_(combs, n = ~n())
  combs <- dplyr::arrange_(combs, ~desc(n))
  res <- unlist(c(unique(combs$V1), utils::tail(combs, 1)[, 2]))
  names(res) <- NULL
  res
}

set_to_formula <- function(x) {
  dep <- x[2]
  ind <- x[1]
  cond <- x[c(-1, -2)]

  stats::formula(paste(dep, paste(c(cond, ind), collapse = '+'), sep = '~'))
}

find_formulas <- function(d, order) {
  s <- ggm::basiSet(d)
  s <- lapply(s, function(x) {
    if (which(order == x[1]) < which(order == x[2])) {
      return(x)
    } else {
      return(c(x[2], x[1], x[-(1:2)]))
    }
  } )
  lapply(s, set_to_formula)
}

C_stat <- function(ps) -2 * sum(log(ps))

C_p <- function(C, k) 1 - stats::pchisq(C, 2 * k)

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

get_p <- function(m) {
  s <- summary(m)$tTable
  s[nrow(s), 'p-value']
}

get_est <- function(m) summary(m)$tTable[-1, 'Value']

get_corStruct <- function(m) attr(m$apVar, "Pars")["corStruct"]