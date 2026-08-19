# The d-separation claims of two models are only comparable when the models
# describe the same set of variables, so a model set has to agree on them.
stop_if_variables_differ <- function(model_set, call = rlang::caller_env()) {
  var_names <- lapply(model_set, colnames)
  if (length(model_set) < 2) return(invisible(model_set))
  if (stats::var(lengths(model_set)) == 0 &&
      all(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) == 0)) {
    return(invisible(model_set))
  }
  shared <- Reduce(intersect, var_names)
  rlang::abort(
    c('All causal models need to include the same variables.',
      x = paste0('Not in every model: ',
                 paste(setdiff(sort(unique(unlist(var_names))), shared), collapse = ', '), '.'),
      i = 'Use the `.common` argument of `define_model_set()` for paths that every model shares.'),
    call = call
  )
}

check_models_data_tree <- function(model_set, data, tree, na.rm) {
  var_names <- lapply(model_set, colnames)
  stop_if_variables_differ(model_set)
  # Check that all variables used in the models are actually present
  used_vars <- unique(unlist(var_names))
  missing_vars <- setdiff(used_vars, colnames(data))
  if (length(missing_vars) > 0) {
    rlang::abort(
      c('Every variable in the causal models has to be a column of `data`.',
        x = paste0('Missing from `data`: ', paste(missing_vars, collapse = ', '), '.')),
      call = rlang::caller_env()
    )
  }
  # Drop all data columns that are not used in the models
  data <- data[, used_vars, drop = FALSE]
  data <- as_binary_factors(data)
  # Check for data columns that are numeric, but only have 0-1. These might be user error, attempting
  # to pass a binary variable as numeric
  binary_numerics <- colnames(data)[sapply(data, function(x) all(x %in% 0:1))]
  if (length(binary_numerics) > 0) {
    rlang::warn(
      c('Some columns contain only 0 and 1, but are numeric, so they are modelled as continuous.',
        '*' = paste0('Columns: ', paste(binary_numerics, collapse = ', '), '.'),
        i = 'Convert them to a factor if they are meant to be binary.')
    )
  }
  # Check tree
  if (inherits(tree, 'multiPhylo')) {
    rlang::abort(
      c('`tree` must be a single phylogeny, but is a `multiPhylo` of several.',
        i = 'Consider using the maximum clade credibility tree, see `?phangorn::mcc`.'),
      call = rlang::caller_env()
    )
  }
  if (!inherits(tree, 'phylo')) {
    rlang::abort(
      c('`tree` must be a phylogeny of class `phylo`.',
        x = paste0('It is of class ', paste(class(tree), collapse = '/'), '.')),
      call = rlang::caller_env()
    )
  }
  # Check NAs and if models and tree line up
  if (anyNA(data)) {
    if (na.rm) {
      NAs <- which(apply(data, 1, anyNA))
      rlang::inform(paste(length(NAs), 'rows were dropped because they contained NA values.'))
      data <- data[-NAs, ]
    } else {
      rlang::abort(
        c('The variables used in the causal models contain NA values.',
          i = 'Set `na.rm = TRUE` to drop the species with missing data, or supply complete data.'),
        call = rlang::caller_env()
      )
    }
  }
  # Match the tree
  unmatched <- setdiff(rownames(data), tree$tip.label)
  if (length(unmatched) > 0) {
    shown <- utils::head(unmatched, 5)
    rlang::abort(
      c('Every species in `data` must be a tip of `tree`, matched exactly by name.',
        x = paste0(length(unmatched), ' of the ', nrow(data),
                   ' rownames of `data` have no matching tip: ',
                   paste(shown, collapse = ', '),
                   if (length(unmatched) > length(shown)) ', ...' else '', '.'),
        i = 'Species names are taken from the rownames of `data`. Tip labels often separate genus and species with an underscore where data uses a space.'),
      call = rlang::caller_env()
    )
  }
  # Prune the tree
  if (length(tree$tip.label) > nrow(data)) {
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(data)))
    rlang::inform('Pruned tree to drop species not included in `data`.')
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
      rlang::abort(
        c(paste0('`', names(data)[i], '` is expected to be binary, but has ', n_levels,
                 ' levels.'),
          x = paste0('Levels: ', paste(levels(data[[i]]), collapse = ', '), '.')),
        call = rlang::caller_env()
      )
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
    rlang::abort(
      c(paste0('Cannot automatically name more than ', length(labels), ' causal models.'),
        x = paste0('The model set has ', n, ' models.'),
        i = 'Name them explicitly, as in `define_model_set(my_name = ...)`.'),
      call = rlang::caller_env()
    )
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
    rlang::abort(
      c('Each causal model needs a unique name.',
        x = paste0('Duplicated: ', paste(unique(nms[duplicated(nms)]), collapse = ', '), '.'),
        i = 'Models are selected by name with `best()` and `choice()`, so names cannot repeat.'),
      call = rlang::caller_env()
    )
  }
  nms
}

# `best()` and `average()` both rank models by CICc, which is meaningless if it
# could not be calculated for some of them.
stop_if_no_CICc <- function(d, fun) {
  if (anyNA(d$CICc)) {
    rlang::abort(
      c(paste0('`', fun, '()` compares the causal models by CICc, so it needs a value for ',
               'every model in the set.'),
        x = paste0('CICc could not be calculated for: ',
                   paste(d$model[is.na(d$CICc)], collapse = ', '), '.')),
      call = rlang::caller_env()
    )
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
    rlang::abort('More than one variable was left out of the consensus ordering. (code 2)',
                 .internal = TRUE)
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
    rlang::abort(
      c('A causal model is fully connected, so it has no independence claims to test.',
        i = 'Phylogenetic path analysis compares models by the independencies they imply, so a model with a path between every pair of variables cannot be evaluated.'),
      call = rlang::caller_env()
    )
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
      rlang::abort('The independence claims came back in an unexpected orientation. (code 1)',
                   .internal = TRUE)
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
    rlang::warn(
      c('Some arguments are not recognized by `phylolm::phylolm()` or `phylolm::phyloglm()`, and are ignored.',
        x = paste0('Ignored: ', paste(unknown, collapse = ', '), '.'))
    )
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
    rlang::abort(
      c('`labels` must be a named vector.',
        i = 'Give it the form `c(old_name = "new label", ...)`.'),
      call = rlang::caller_env()
    )
  }
  if (length(setdiff(l$name, names(labels))) > 0) {
    rlang::abort(
      c('`labels` must name every variable in the causal model.',
        x = paste0('Missing from `labels`: ',
                   paste(setdiff(l$name, names(labels)), collapse = ', '), '.')),
      call = rlang::caller_env()
    )
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
  rlang::warn(
    c('Cannot determine which variables of this model are binary, so the scale of its coefficients is unknown.',
      ' ' = 'Paths into a binary variable are log odds ratios, while paths into a continuous variable are standardized regression coefficients.',
      i = 'This model was fitted by a version of phylopath older than 1.4.0, which did not record it. Refit it to label the coefficients, and to stop receiving this warning.')
  )
}

# The variables of a fitted model, ordered from upstream to downstream.
causal_order <- function(coef) {
  adjacency <- (coef != 0) * 1
  if (!ggm::isAcyclic(adjacency)) {
    # The variables caught in a cycle are those in a strongly connected component
    # of more than one variable. That covers a pair of paths pointing both ways at
    # each other as well as a longer loop, both of which averaging can produce.
    comp <- igraph::components(igraph::graph_from_adjacency_matrix(adjacency),
                               mode = 'strong')
    in_cycle <- names(comp$membership)[comp$membership %in% which(comp$csize > 1)]
    rlang::abort(
      c('`order_by = "causal"` follows the paths of a causal model from upstream to downstream, which is only possible when the model is acyclic.',
        x = paste0('The paths between ', paste(in_cycle, collapse = ', '), ' form a cycle.'),
        ' ' = 'Averaging can produce a cyclical model, when the direction of a path is not resolved and different models in the set point it different ways.',
        i = 'Use `order_by = "default"` or `"strength"` instead.'),
      call = rlang::caller_env()
    )
  }
  colnames(ggm::topSort(adjacency))
}

# Get the variables that take part in at least one path.
used_vars <- function(coef) {
  present <- coef != 0
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
    rlang::abort(
      c('The fitted models disagree about which variables are binary.',
        x = paste0('Modelled both ways: ', paste(disagree, collapse = ', '), '.'),
        i = 'Coefficients of paths into binary and continuous variables are not on the same scale, so they cannot be averaged together.'),
      call = rlang::caller_env()
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
    rlang::abort('`x`, `se` and `weight` must be numeric vectors.',
                 call = rlang::caller_env())
  n <- length(x)
  if (length(weight) != n || length(se) != n) {
    rlang::abort(
      c('`x`, `se` and `weight` must all be the same length.',
        x = paste0('They are ', length(x), ', ', length(se), ' and ', length(weight),
                   ' long.')),
      call = rlang::caller_env()
    )
  }
  weight[is.na(weight)] <- 0
  wx <- stats::weighted.mean(x, weight, na.rm = TRUE)

  x.sqdiff <- (x - wx) ^ 2
  xvar <- se ^ 2

  se <- sqrt(stats::weighted.mean(xvar + x.sqdiff, weight, na.rm = TRUE))
  ci <- stats::qnorm(0.975, lower.tail = TRUE) * se

  return(c(Coefficient = wx, SE = se, `Lower CI` = wx - ci, `Upper CI` = wx + ci))
}

# Paths named as "from -> to", so that coef() and confint() line up.
path_names <- function(paths) {
  paste(paths$from, '->', paths$to)
}
