check_models_data_tree <- function(models, data, tree, na.rm) {
  var_names <- lapply(models, colnames)
  if (length(models) > 1 &
      (stats::var(lengths(models)) != 0 |
       any(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) != 0))) {
    stop('All causal models need to include the same variables. Combined, your
       models include the following variables:\n',
         paste(sort(unique(unlist(var_names))), collapse = '\n'),
         call. = FALSE)
  }
  data <- data[, unique(unlist(var_names))]
  # Check NAs and if models and tree line up
  if ('tbl_df' %in% class(data)) {
    data <- as.data.frame(data)
  }
  if (anyNA(data)) {
    if (na.rm) {
      NAs <- which(apply(data, 1, anyNA))
      message(length(NAs), ' rows were dropped because they contained NA values.')
      data <- data[-NAs, ]
    } else {
      stop('NA values were found in the variables of interest.', call. = FALSE)
    }
  }
  if (length(setdiff(rownames(data), tree$tip.label)) > 0) {
    stop('Make sure that species in your data have rownames that are exactly matched by name with tips in the tree.')
  }
  if (length(tree$tip.label) > nrow(data)) {
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(data)))
    message('Pruned tree to drop species not included in dat.')
  }
  if (is.null(names(models))) {
    names(models) <- LETTERS[1:length(models)]
  }
  return(list(models = models, data = data, tree = tree))
}

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
    m <- purrr::safely(function(.x, ...) nlme::gls(..., correlation = cor_fun(par, .x)))(tree, ...)
    if (is.null(m$error)) break
  }
  return(m)
}

get_p <- function(m) {
  s <- stats::coef(summary(m))
  s[nrow(s), ncol(s)]
}

get_p_binary <- function(m) {
  s <- m$B.pvalue[, 1]
  utils::tail(s, 1)
}

get_est <- function(m) {
  s <- stats::coef(summary(m))
  s[-1, 1]
}

get_est_binary <- function(m) {
  m$B[-1, 1]
}

get_se <- function(m) {
  s <- stats::coef(summary(m))
  s[-1, 2]
}

get_se_binary <- function(m) {
  m$B.se[-1, 1]
}

get_lower <- function(m) nlme::intervals(m)$coef[-1, 'lower']

get_upper <- function(m) nlme::intervals(m)$coef[-1, 'upper']

get_corStruct <- function(m) m$modelStruct

adjust_layout <- function(l, rotation, flip_x, flip_y) {
  rotation <- rotation * (2 * pi / 360)
  R <- matrix(c(cos(rotation), sin(rotation), -sin(rotation), cos(rotation)), nrow = 2)
  l[c('x', 'y')] <- as.matrix(l[c('x', 'y')]) %*% R
  if (flip_x) {
    l$x <- -l$x
  }
  if (flip_y) {
    l$y <- -l$y
  }
  return(l)
}

combine_with_labels <- function(l, labels) {
  if (is.null(labels)) {
    return(l)
  }
  if (is.null(names(labels))) {
    stop('labels must be a named vector.', call. = FALSE)
  }
  if (length(setdiff(l$name, names(labels))) > 0) {
    stop('Some nodes are missing from labels.', call. = FALSE)
  }
  l$name <- factor(l$name, names(labels), labels)
  class(l) <- c("layout_igraph", "layout_ggraph", "data.frame")
  return(l)
}

#' @importFrom ggraph guide_train.edge_colourbar
#' @export
ggraph::guide_train.edge_colourbar
