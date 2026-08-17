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
  # Check that all variables used in the models are actually present
  used_vars <- unique(unlist(var_names))
  missing_vars <- setdiff(used_vars, colnames(data))
  if (length(missing_vars) > 0) {
    stop('These variables are used in the causal models, but are not columns in `data`: ',
         paste(missing_vars, collapse = ', '), '.', call. = FALSE)
  }
  # Drop all data columns that are not used in the models
  data <- data[, used_vars, drop = FALSE]
  data <- as_binary_factors(data)
  # Check for data columns that are numeric, but only have 0-1. These might be user error, attempting
  # to pass a binary variable as numeric
  binary_numerics <- sapply(data, function(x) all(x %in% 0:1))
  for (n in colnames(data)[binary_numerics]) {
    warning('Column ', n, ' appears to have binary data, but was not recognized as binary. If it should be treated as binary, convert it to a factor first.')
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
    message('Pruned tree to drop species not included in `data`.')
  }
  # Add names to any models that don't have them
  names(model_set) <- name_model_set(names(model_set), length(model_set))
  return(list(model_set = model_set, data = data, tree = tree))
}

# Binary variables can be passed as either characters or factors, but only
# factors are used internally. Characters are therefore converted here,
# with some input checks.
as_binary_factors <- function(data) {
  char_cols <- sapply(data, is.character)
  data[char_cols] <- lapply(data[char_cols], as.factor)
  # Check whether all factors have exactly two levels:
  f_cols <- which(sapply(data, is.factor))
  for (i in f_cols) {
    n_levels <- length(levels(data[[i]]))
    if (n_levels != 2) {
      stop("Variable '", names(data)[i], "' is expected to be binary, but has ", n_levels,
           " levels.", call. = FALSE)
    }
  }
  data
}

# Labels for automatically named models: A, B, ... Z, then AA, AB, ... ZZ.
model_set_labels <- function(n) {
  labels <- LETTERS
  if (n > length(labels)) {
    labels <- c(labels, paste0(rep(LETTERS, each = length(LETTERS)), LETTERS))
  }
  if (n > length(labels)) {
    stop('Cannot automatically name more than ', length(labels),
         ' causal models. Please name your models explicitly.', call. = FALSE)
  }
  labels[seq_len(n)]
}

# Fill in names for any models the user left unnamed. Models are looked up by
# name downstream (by `best()` and `choice()`), so every model needs a name and
# all names have to be unique.
name_model_set <- function(nms, n) {
  if (is.null(nms)) nms <- rep('', n)
  blank <- is.na(nms) | nms == ''
  if (any(blank)) {
    nms[blank] <- model_set_labels(n)[blank]
  }
  if (anyDuplicated(nms)) {
    stop('Each causal model needs a unique name, but the following are duplicated: ',
         paste(unique(nms[duplicated(nms)]), collapse = ', '), '.', call. = FALSE)
  }
  nms
}

# `best()` and `average()` both rank models by CICc, which is meaningless if it
# could not be calculated for some of them.
stop_if_no_CICc <- function(d, fun) {
  if (anyNA(d$CICc)) {
    stop('CICc could not be calculated for the following causal model(s): ',
         paste(d$model[is.na(d$CICc)], collapse = ', '), '.\n',
         '  `', fun, '()` compares models by CICc, so it needs a value for every ',
         'model in the set. Either use a smaller set of models, or collect data ',
         'for more species (CICc requires n > q + 1).', call. = FALSE)
  }
  invisible(d)
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
  # find the unique variables in the models
  vars <- lapply(model_set, colnames)
  # tally each possible combination of variables
  combs <- as.data.frame(t(utils::combn(vars[[1]], 2)), stringsAsFactors = FALSE)
  names(combs) <- c('node1', 'node2')
  # count how often var1 is ordered above var2
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

  # Now we order the nodes by how many nodes they are above (i.e. how often are they node 1?),
  # this should go from n:1
  combs$n <- table(combs$node1)[combs$node1]
  combs <- combs[order(-combs$n), ]
  # now we find all the node1 variables, in order of n
  res <- unique(combs$node1)
  # but we might be missing one last variable, if it was *always* node 2
  missing_vars <- setdiff(vars[[1]], res)
  # If so, add at the end of the order
  if (length(missing_vars) == 1) {
    res <- c(res, missing_vars)
  }
  # we should never have more than 1 missing var, but check for this just in case:
  if (length(missing_vars) > 1) {
    stop('If you get this error, please contact the maintainer. (code 2)')
  }

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
      stop('If you get this error, please contact the maintainer. (code 1)')
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

# Computed as the upper tail directly: 1 - pchisq() underflows to exactly 0
# once the tail is smaller than about 1e-16, and 0 is never a correct p-value.
C_p <- function(C, k) stats::pchisq(C, 2 * k, lower.tail = FALSE)

CICc <- function(C, q, n) C + 2 * q * (n / (n - 1 - q))

l <- function(dCICc) exp(-0.5 * dCICc)

w <- function(l) l / sum(l)

# Drop any arguments in ... that neither phylolm nor phyloglm accepts, and generate
# a warning when they are passed.
check_dots <- function(dots) {
  known <- union(names(formals(phylolm::phyloglm)), names(formals(phylolm::phylolm)))
  unknown <- setdiff(names(dots), known)
  if (length(unknown) > 0) {
    warning('The following arguments are not recognized by `phylolm::phylolm()` or ',
            '`phylolm::phyloglm()`, and are ignored: ', paste(unknown, collapse = ', '),
            call. = FALSE)
  }
  dots[!(names(dots) %in% unknown)]
}

phylo_g_lm <- function(formula, data, tree, model, method, boot = 0, dots = NULL) {
  # we capture the dots, because we need to match the names to either phylolm or phylolm
  dots_glm <- dots[names(dots) %in% names(formals(phylolm::phyloglm))]
  dots_lm <- dots[names(dots) %in% names(formals(phylolm::phylolm))]
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
              dots_lm)
  }
  res <- do.call(quiet_safely(fun), args)
  # Remove the call, since quiet_safely messes it up and it's annoying in printing
  res$result$call <- NULL

  return(res)
}

get_p <- function(m) {
  s <- switch(
    class(m),
    phylolm = stats::coef(phylolm::summary.phylolm(m)),
    phyloglm = stats::coef(phylolm::summary.phyloglm(m))
  )
  p <- s[nrow(s), 'p.value']
  # Floor only at the smallest representable number, so that C = -2 * sum(log(p))
  # cannot become infinite, without discarding p-values that are genuinely small:
  # phylolm computes these as 2 * pt(-abs(t), df), which stays accurate well
  # below machine epsilon.
  if (p < .Machine$double.xmin) p <- .Machine$double.xmin
  return(p)
}

get_est <- function(m) {
  stats::coef(m)[-1]
}

get_se <- function(m) {
  s <- switch(
    class(m),
    phylolm = stats::coef(phylolm::summary.phylolm(m)),
    phyloglm = stats::coef(phylolm::summary.phyloglm(m))
  )
  s[-1, 'StdErr']
}

get_lower <- function(m) {
  s <- switch(
    class(m),
    phylolm = stats::coef(phylolm::summary.phylolm(m)),
    phyloglm = stats::coef(phylolm::summary.phyloglm(m))
  )
  if ('lowerbootCI' %in% colnames(s)) {
    r <- s[-1, 'lowerbootCI']
  } else {
    r <- NA
  }
  return(r)
}

get_upper <- function(m) {
  s <- switch(
    class(m),
    phylolm = stats::coef(phylolm::summary.phylolm(m)),
    phyloglm = stats::coef(phylolm::summary.phyloglm(m))
  )
  if ('upperbootCI' %in% colnames(s)) {
    r <- s[-1, 'upperbootCI']
  } else {
    r <- NA
  }
  return(r)
}

get_phylo_param <- function(m) {
  r <- m$optpar
  # phyloglm reports its phylogenetic parameter as `alpha`, not as `optpar`.
  if (is.null(r)) r <- m$alpha
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
  incoming_class <- class(l)
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
  class(l) <- incoming_class
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

# Warning thrown whenever a function has no information on compatiable scales.
warn_unknown_scale <- function() {
  warning('Cannot determine which variables of this model are binary, so the scale of ',
          'its coefficients is unknown: paths into a binary variable are log odds ',
          'ratios, while paths into a continuous variable are standardized regression ',
          'coefficients. This model was fitted by a version of phylopath older than ',
          '1.4.0, which did not record it. Refit the model to label the coefficients, ',
          'and stop receiving these warnings.',
          call. = FALSE)
}

# Get the variables that take part in at least one path.
used_vars <- function(coef) {
  present <- abs(coef) > .Machine$double.eps
  union(rownames(coef)[rowSums(present) > 0], colnames(coef)[colSums(present) > 0])
}

# Is this a model whose coefficients cannot be compared with each other? That is
# the case when the variables it connects are not all binary or all continuous.
# NA when the variable types are not recorded.
is_mixed <- function(fitted_DAG) {
  binary <- fitted_DAG$binary
  coef <- fitted_DAG$coef
  if (is.null(binary) || !all(colnames(coef) %in% names(binary))) return(NA)

  return(length(unique(binary[used_vars(coef)])) > 1)
}

# The scale that every coefficient of this model is on, or NULL when the model
# mixes types and there is no single answer. Warns when the variable types were
# not recorded, in which case nothing can be said about them at all. Warning here
# rather than at each call site means every function that reports coefficients
# does so consistently.
coef_scale_name <- function(fitted_DAG) {
  mixed <- is_mixed(fitted_DAG)
  if (is.na(mixed)) {
    warn_unknown_scale()
    return(NULL)
  }
  if (mixed) return(NULL)
  if (any(fitted_DAG$binary[used_vars(fitted_DAG$coef)])) {
    return('log odds ratio')
  } else {
    return('standardized regression coefficient')
  }
}

# Warn against the interpretation of coefficients in mixed DAGs, one per session.
warn_if_mixed <- function(fitted_DAG) {
  if (!isTRUE(is_mixed(fitted_DAG))) return(invisible(FALSE))
    rlang::warn(
      c('This causal model contains both binary and continuous variables, so its path coefficients cannot be compared with one another.',
        ' ' = 'Paths into a binary variable are log odds ratios, while paths into a continuous variable are standardized regression coefficients. On top of that, a binary predictor contributes a contrast between its two levels, which is always more than two standard deviations, so its coefficients are inflated relative to those of continuous predictors.',
        ' ' = 'Comparison of the causal models themselves is not affected: CICc, the model weights and the ranking of models do not depend on the scale of the coefficients.',
        i = 'See `?est_DAG` for what each coefficient means.'),
      .frequency = 'once', .frequency_id = 'phylopath_mixed_type'
    )
  invisible(TRUE)
}

# The `binary` flags shared by a list of fitted_DAGs. Averaging coefficients of
# different types would be meaningless, so disagreement generates an error.
# Returns NULL if the flags are unavailable, which happens for fitted_DAG objects
# made before v1.4.0.
combine_binary <- function(fitted_DAGs, ord) {
  binary <- lapply(fitted_DAGs, function(l) l$binary[ord])
  if (
    any(vapply(binary, function(b) is.null(b) || anyNA(b), logical(1)))
  ) {
    warn_unknown_scale()
    return(NULL)
  }
  if (length(unique(binary)) > 1) {
    disagree <- names(which(apply(simplify2array(binary), 1, function(x) {
      length(unique(x)) > 1
    })))
    stop(
      'The fitted models disagree about which variables are binary: ',
      paste(disagree, collapse = ', '),
      '.\n',
      '  Coefficients of paths into binary and continuous variables are not on ',
      'the same scale, so they cannot be averaged together.',
      call. = FALSE
    )
  }
  binary[[1]]
}

# The polite restatement, for a plot caption or a note below a printed table.
mixed_note <- function() {
  paste('This model mixes binary and continuous variables, so its coefficients',
        'cannot be compared\nwith one another. See ?est_DAG for what each one means.')
}

par_avg <- function(x, se, weight) {
  # Derived from original MuMIn::par.avg function, written by Kamil Bartoń.

  if (!(is.numeric(x) && is.numeric(se) && is.numeric(weight)))
    stop("'x', 'se' and 'weight' must be numeric vectors")
  n <- length(x)
  if (length(weight) != n || length(se) != n) {
    stop("'x', 'se' and 'weight' are not of the same length")
  }
  weight[is.na(weight)] <- 0
  wx <- stats::weighted.mean(x, weight, na.rm = TRUE)

  x.sqdiff <- (x - wx) ^ 2
  xvar <- se ^ 2

  se <- sqrt(stats::weighted.mean(xvar + x.sqdiff, weight, na.rm = TRUE))
  ci <- stats::qnorm(0.975, lower.tail = TRUE) * se

  return(c(Coefficient = wx, SE = se, `Lower CI` = wx - ci, `Upper CI` = wx + ci))
}
