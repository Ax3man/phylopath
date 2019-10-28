check_models_data_tree <- function(model_set, data, tree, na.rm) {
  var_names <- lapply(model_set, colnames)
  # Check whether all causal models have the same variables.
  if (length(model_set) > 1 &
      (stats::var(lengths(model_set)) != 0 |
       any(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) != 0))) {
    stop('All causal models need to include the same variables. Your
       model set includes the following variables:\n',
         paste(sort(unique(unlist(var_names))), collapse = '\n'),
         call. = FALSE)
  }
  data <- data[, unique(unlist(var_names))]
  # We force all character columns to factors
  char_cols <- sapply(data, is.character)
  data[char_cols] <- lapply(data[char_cols], as.factor)
  # Check whether all factors have exactly two levels:
  f_cols <- which(sapply(data, is.factor))
  for (i in f_cols) {
    n_levels <- length(levels(data[[i]]))
    if (n_levels != 2) {
      stop("Variable '", names(data)[i], "' is expected to binary, but has ", n_levels, " levels.",
           .call = FALSE)
    }
  }
  # Check tree
  if (inherits(tree, 'multiPhylo')) {
    stop('You are passing several trees (in a `multiPhylo` object). Please only pass one `phylo` object.')
  }
  if (!inherits(tree, 'phylo')) {
    stop('The tree needs to be of class `phylo`.')
  }
  # Check NAs and if models and tree line up
  if (anyNA(data)) {
    if (na.rm) {
      NAs <- which(apply(data, 1, anyNA))
      message(length(NAs), ' rows were dropped because they contained NA values.')
      data <- data[-NAs, ]
    } else {
      stop('NA values were found in the variables of interest.', call. = FALSE)
    }
  }
  # Match the tree
  if (length(setdiff(rownames(data), tree$tip.label)) > 0) {
    stop('Make sure that species in your data have rownames that are exactly matched by name with tips in the tree.')
  }
  # Prune the tree
  if (length(tree$tip.label) > nrow(data)) {
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(data)))
    message('Pruned tree to drop species not included in dat.')
  }
  # Add names to the models, if they don't have them
  if (is.null(names(model_set))) {
    names(model_set) <- LETTERS[1:length(model_set)]
  }
  return(list(model_set = model_set, data = data, tree = tree))
}

find_consensus_order <- function(model_set) {
  # If the fully combined model is acyclic, then we use that.
  model_set_same_order <- lapply(model_set, function(x) {
    x[rownames(model_set[[1]]), colnames(model_set[[1]])]
  } )
  full_model <- sign(Reduce('+', model_set_same_order))
  if (ggm::isAcyclic(full_model)) {
    return(rownames(ggm::topSort(full_model)))
  }
  # Otherwise we find the most common orderings and use those.
  # Make sure all models are ordered:
  model_set <- lapply(model_set, ggm::topSort)
  vars <- lapply(model_set, colnames)
  combs <- as.data.frame(t(utils::combn(vars[[1]], 2)), stringsAsFactors = FALSE)
  names(combs) <- c('node1', 'node2')
  combs$count <- 0
  for (i in seq_along(vars)) {
    v <- apply(combs, 1, function(x) {
      which(vars[[i]] == x[1]) < which(vars[[i]] == x[2])
    } )
    combs$count <- combs$count + v
  }

  # If node1 is commonly ordered above node2, leave as is, otherwise swap them around
  tmp <- combs$node1
  combs$node1 <- ifelse(combs$count > 0.5 * length(model_set), combs$node1, combs$node2)
  combs$node2 <- ifelse(combs$count > 0.5 * length(model_set), combs$node2, tmp)

  # Now we order the nodes by how many nodes they are above, this should go from n:1
  combs$n <- table(combs$node1)[combs$node1]
  combs <- combs[order(-combs$n), ]
  res <- unlist(c(unique(combs$node1), utils::tail(combs, 1)[, 'node2']))
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
  if (is.null(s)) {
    stop('One or some of your models are fully connected, and cannot be tested.')
  }
  s <- lapply(s, function(x) {
    # define whether there are existing paths between the two nodes in both directions.
    path1 <- !is.null(ggm::findPath(d, which(rownames(d) == x[1]), which(rownames(d) == x[2])))
    path2 <- !is.null(ggm::findPath(d, which(rownames(d) == x[2]), which(rownames(d) == x[1])))
    if (path1 & !path2) {
      # the first vertex is upstream, so we do not re-order
      return(x)
    }
    if ((path2 & !path1) | (path1 & path2)) {
      # these conditions should not occur, the first means basiSet is returning the wrong order,
      # the second should only occur if there are cycles.
      stop('If you get this error, please contact the maintainer.')
    }
    if (!path1 & !path2) {
      # check whether the order is according to `order`
      if (which(order == x[1]) < which(order == x[2])) {
        return(x)
      } else {
        return(c(x[2], x[1], x[-(1:2)]))
      }
    }
  } )
  lapply(s, set_to_formula)
}

C_stat <- function(ps) -2 * sum(log(ps))

C_p <- function(C, k) 1 - stats::pchisq(C, 2 * k)

CICc <- function(C, q, n) C + 2 * q * (n / (n - 1 - q))

l <- function(dCICc) exp(-0.5 * dCICc)

w <- function(l) l / sum(l)

phylo_g_lm <- function(formula, data, tree, model, method, boot = 0, ...) {
  # we capture the dots, because we need to match the names to either phylolm or phylolm
  dots <- list(...)
  dots_glm <- dots[names(dots) %in% names(formals(phylolm::phyloglm))]
  dots_lm <- dots[names(dots) %in% names(formals(phylolm::phylolm))]
  if (length(intersect(names(dots_glm), names(dots_lm))) != length(dots)) {
    warning("Some arguments in ... are not recognized.", call. = FALSE)
  }
  # we capture the first argument in the formula, to check whether it is binary
  x_var <- data[[all.vars(formula)[1]]]
  if (is.factor(x_var)) {
    # phyloglm need binary variables as 0,1 but I use factors
    data[all.vars(formula)[1]] <- as.numeric(x_var) - 1
    fun <- phylolm::phyloglm
    args <- c(list(formula = formula, data = data, phy = tree, method = method, boot = boot),
              dots_glm)
  } else {
    fun <- phylolm::phylolm
    args <- c(list(formula = formula, data = data, phy = tree, model = model, boot = boot),
              dots_glm)
  }
  res <- do.call(quiet_safely(fun), args)
  # Remove the call, since quiet_safely messes it up and it's annoying in printing
  res$result$call <- NULL

  return(res)
}

get_p <- function(m) {
  s <- stats::coef(summary(m))
  p <- s[nrow(s), 'p.value']
  if (p < .Machine$double.eps) p <- .Machine$double.eps
  return(p)
}

get_est <- function(m) {
  stats::coef(m)[-1]
}

get_se <- function(m) {
  stats::coef(summary(m))[-1, 'StdErr']
}

get_lower <- function(m) {
  s <- stats::coef(summary(m))
  if ('lowerbootCI' %in% colnames(s)) {
    r <- s[-1, 'lowerbootCI']
  } else {
    r <- NA
  }
  return(r)
}

get_upper <- function(m) {
  s <- stats::coef(summary(m))
  if ('upperbootCI' %in% colnames(s)) {
    r <- s[-1, 'upperbootCI']
  } else {
    r <- NA
  }
  return(r)
}

get_phylo_param <- function(m) {
  r <- m$optpar
  if (is.null(r)) r <- NA
  return(r)
}

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

quiet_safely <- function(.f) {
  capture_all <- function(expr)
  {
    warn_vec <- NULL
    w.handler <- function(w){ # warning handler
      warn_vec <<- c(warn_vec, w$message)
      invokeRestart("muffleWarning")
    }
    r <- list(result = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
              warning = warn_vec)
    if (inherits(r$result, 'error')) {
      return(list(result = NULL, error = r$result$message, warning = r$warning))
    } else {
    return(list(result = r$result, error = NULL, warning = r$warning))
    }
  }
  function(...) capture_all(.f(...))
}

combine_dots <- function(old_dots, ...) {
  new_dots <- list(...)
  c(new_dots, old_dots[!(names(old_dots) %in% names(new_dots))])
}
