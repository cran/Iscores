#' Combine two projection forests
#'
#' @param mod1 A fitted forest object.
#' @param mod2 A fitted forest object.
#'
#' @return A forest object containing trees from both input forests.
#' @keywords internal

combine2Forests <- function(mod1, mod2) {
  res <- mod1

  res$num.trees <- mod1$num.trees + mod2$num.trees
  res$inbag.counts <- c(mod1$inbag.counts, mod2$inbag.counts)
  res$forest$child.nodeIDs <- c(mod1$forest$child.nodeIDs, mod2$forest$child.nodeIDs)

  for (i in seq_len(mod1$num.trees)) {
    res$forest$split.varIDs[[i]] <- mod1$var[mod1$forest$split.varIDs[[i]] + 1]
    res$forest$split.varIDs[[i]][mod1$forest$child.nodeIDs[[i]][[1]] == 0] <- 0
  }

  for (i in seq_len(mod2$num.trees)) {
    mod2$forest$split.varIDs[[i]] <- mod2$var[mod2$forest$split.varIDs[[i]] + 1]
    mod2$forest$split.varIDs[[i]][mod2$forest$child.nodeIDs[[i]][[1]] == 0] <- 0
  }

  res$forest$split.varIDs <- c(res$forest$split.varIDs, mod2$forest$split.varIDs)
  res$forest$split.values <- c(mod1$forest$split.values, mod2$forest$split.values)

  if (!is.null(mod1$forest$terminal.class.counts)) {
    res$forest$terminal.class.counts <- c(
      mod1$forest$terminal.class.counts,
      mod2$forest$terminal.class.counts
    )
  }

  res$forest$num.trees <- mod1$forest$num.trees + mod2$forest$num.trees
  res$call$num.trees <- res$num.trees
  res$num.independent.variables <- length(res$full.vars)
  res$forest$independent.variable.names <- res$full.vars
  res$forest$is.ordered <- rep(TRUE, length(res$full.vars))
  res$var <- 0:(length(res$full.vars) - 1)

  res
}


#' Combine a list of forests
#'
#' @param list.rf A list of fitted forest objects.
#'
#' @return A single forest object obtained by combining all forests in \code{list.rf}.
#' @keywords internal

combineForests <- function(list.rf) {
  Reduce(combine2Forests, list.rf)
}


#' Computation of the density ratio score
#'
#' @importFrom kernlab kernelMatrix
#' @importFrom kernlab rbfdot
#'
#' @description Computes the density ratio score using a random forest model
#' based on random projections.
#'
#' @param X A numeric matrix of observed data that may contain missing values
#' denoted by \code{NA}.
#' @param X_imp A numeric matrix of imputed values with the same dimensions as
#' \code{X}.
#' @param pattern A vector or pattern indicating the missingness structure.
#' @param n_proj An integer specifying the number of random projections.
#' @param n_trees_per_proj An integer specifying the number of trees grown per
#' projection.
#' @param projection_function A function that generates user-defined projections.
#' @param min_node_size An integer specifying the minimum number of observations
#' in a terminal node (leaf) of each tree.
#' @param normal_proj Logical. If \code{TRUE}, sampling is performed from both
#' missing (NA) and observed values. If \code{FALSE}, sampling is performed only
#' from missing (NA) values.
#'
#' @return An object representing a fitted random forest model based on random
#' projections.
#'
#' @details The method builds multiple random forests on projected versions of
#' the data to estimate the density ratio between observed and imputed
#' distributions.
#'
#' @keywords internal

densityRatioScore <- function(X,
                              X_imp,
                              pattern = NULL,
                              n_proj = 10,
                              n_trees_per_proj = 1,
                              projection_function = NULL,
                              min_node_size = 1,
                              normal_proj = TRUE) {
  M <- is.na(as.matrix(X))

  if (!is.null(pattern)) {
    ids_x_na <- which(as.logical(pattern))
  } else {
    ids_x_na <- 0
  }

  list_rf <- lapply(seq_len(n_proj), function(i) {
    vars <- sample_vars_proj(
      ids_x_na = ids_x_na,
      X = X,
      projection_function = projection_function,
      normal_proj = normal_proj
    )

    dim_proj <- length(vars)

    X_proj_complete <- as.matrix(stats::na.omit(X[, vars, drop = FALSE]))
    ids_with_missing <- which(rowSums(M[, vars, drop = FALSE]) != 0)

    if (nrow(X_proj_complete) <= 2 || length(ids_with_missing) == 0) {
      return(NA)
    }

    if (normal_proj) {
      patternxA <- matrix(pattern[vars], ncol = length(vars))
      M_A <- M[ids_with_missing, vars, drop = FALSE]
      kern <- kernlab::rbfdot(sigma = 0.25)
      B <- kernlab::kernelMatrix(kern, x = patternxA, y = M_A)
      drawA <- which(B == 1)
      Y_proj <- X_imp[ids_with_missing, vars, drop = FALSE][drawA, , drop = FALSE]
    } else {
      Y_proj <- X_imp[ids_with_missing, vars, drop = FALSE]
      drawA <- integer(0)
    }

    cl_bl_output <- class.balancing(
      X_proj_complete = X_proj_complete,
      Y.proj = Y_proj,
      drawA = drawA,
      X_imp = X_imp,
      ids.with.missing = ids_with_missing,
      vars = vars
    )

    X_proj_complete <- cl_bl_output$X_proj_complete
    Y_proj <- cl_bl_output$Y.proj

    if (nrow(Y_proj) == 0) {
      return(NA)
    }

    colnames(Y_proj) <- NULL
    Y_proj <- as.matrix(Y_proj)
    colnames(X_proj_complete) <- NULL

    d <- data.frame(
      class = c(
        rep(1, each = nrow(X_proj_complete)),
        rep(0, each = nrow(Y_proj))
      ),
      X = rbind(X_proj_complete, Y_proj)
    )

    obj <- tryCatch(
      {
        ranger::ranger(
          probability = TRUE,
          formula = class ~ .,
          data = d,
          num.trees = n_trees_per_proj,
          mtry = dim_proj,
          keep.inbag = TRUE,
          min.node.size = min_node_size
        )
      },
      error = function(e) NA
    )

    if (inherits(obj, "try-error")) {
      return(NA)
    }

    if (any(is.na(obj))) {
      warning("Forest for a projection was NA, will redo")
    } else {
      obj$var <- vars - 1
      obj$full.vars <- paste("X.", 1:ncol(X), sep = "")
      return(obj)
    }
  })

  list_rf <- list_rf[lengths(list_rf) != 1]

  list(list.rf = combineForests(list.rf = list_rf))
}

#' Compute the density ratio score
#'
#' @importFrom stats predict
#'
#' @param object a crf object.
#' @param Z a matrix of candidate points.
#' @param n_proj an integer specifying the number of projections.
#' @param n_trees_per_proj an integer, the number of trees per projection.
#'
#' @return a numeric value, the DR I-Score.
#' @keywords internal

compute_drScore <- function(object, Z = Z, n_trees_per_proj, n_proj) {

  preds_all_f_h <- predict(
    object$list.rf,
    data.frame(X = Z),
    predict.all = TRUE
  )$predictions

  if (length(dim(preds_all_f_h)) == 2) {

    preds_all_f_h <- apply(preds_all_f_h, 1, mean)
    p_f_h <- matrix(preds_all_f_h, nrow = 1)
    dr_f_h <- matrix(apply(p_f_h, 2,
                           function(p) truncProb(p) / (1 - truncProb(p))),
                     ncol = 1)

    dr_f_h[!is.finite(dr_f_h)] <- 0
    kl_f_h <- log(dr_f_h)
    kl_f_h <- colMeans(kl_f_h)
    scoredr <- kl_f_h
  } else {

    if (!n_trees_per_proj > 1) {

      preds_all_f_h <- t(
        apply(preds_all_f_h, 1, cumsum)
      )[ , seq(n_trees_per_proj, n_trees_per_proj * n_proj, n_trees_per_proj),
         drop = FALSE]

      if (ncol(preds_all_f_h) <= 2) {

        preds_all_f_h <- cbind(
          preds_all_f_h[, 1],
          matrix(apply(preds_all_f_h, 1, diff), ncol = 1) ) / n_trees_per_proj
      } else {

        preds_all_f_h <- cbind(
          preds_all_f_h[, 1], t(apply(preds_all_f_h, 1, diff))
        ) / n_trees_per_proj
      }
    }

    if (dim(preds_all_f_h)[1] == 1) {

      p_f_h <- matrix(preds_all_f_h[1, 1, ], nrow = 1)
      dr_f_h <- matrix(
        apply(p_f_h, 2, function(p) truncProb(p) / (1 - truncProb(p))),
        ncol = 1
      )
    } else {

      p_f_h <- preds_all_f_h[, 1, ]
      dr_f_h <- t(apply(p_f_h, 2, function(p) truncProb(p) / (1 - truncProb(p))))
    }

    dr_f_h[!is.finite(dr_f_h)] <- 0
    kl_f_h <- log(dr_f_h)
    kl_f_h <- colMeans(kl_f_h)
    scoredr <- kl_f_h
  }

  scoredr
}

#' Truncation of probability
#'
#' @param p a numeric value between 0 and 1 to be truncated
#'
#' @return a numeric value, the truncated probability.
#'
#' @keywords internal

truncProb <- function(p) {
  pmin(pmax(p, 10^-9), 1 - 10^-9)
}

#' Balancing of Classes
#'
#' @param X_proj_complete matrix with complete projected observations.
#' @param Y.proj matrix with projected imputed observations.
#' @param drawA vector of indices corresponding to current missingness pattern.
#' @param X_imp matrix of full imputed observations.
#' @param ids.with.missing vector of indices of observations with missing values.
#' @param vars vectors of variables in projection.
#'
#' @return a list of new X_proj_complete and Y.proj.
#' @keywords internal

class.balancing <- function(X_proj_complete,
                            Y.proj,
                            drawA,
                            X_imp,
                            ids.with.missing,
                            vars) {
  if (nrow(Y.proj) >= nrow(X_proj_complete)) {
    X_proj_complete <- X_proj_complete[
      sample(1:nrow(X_proj_complete), size = nrow(Y.proj), replace = TRUE),
    ]
  } else {
    X_imp0 <- X_imp[ids.with.missing, vars, drop = FALSE][-drawA, , drop = FALSE]

    if ((nrow(Y.proj) < nrow(X_proj_complete) * 0.75) && nrow(X_imp0) > 0) {
      ind <- sample(
        1:nrow(X_imp0),
        size = (nrow(X_proj_complete) - nrow(Y.proj)),
        replace = TRUE
      )

      Y.proj <- rbind(Y.proj, X_imp0[ind, ])
    } else {
      Y.proj <- Y.proj[
        sample(1:nrow(Y.proj), size = nrow(X_proj_complete), replace = TRUE),
        ,
        drop = FALSE
      ]
    }
  }

  return(list(X_proj_complete = X_proj_complete, Y.proj = Y.proj))
}

#' Sampling of Projections
#'
#' @param ids_x_na a vector of indices corresponding to NA in the given
#' missingness pattern.
#' @param X a matrix of the observed data containing missing values.
#' @param projection_function a function providing the user-specific projections.
#' @param normal_proj a boolean, if TRUE, sample from the NA of the pattern and
#' additionally from the non-NA. If FALSE, sample only from the NA of the
#' pattern.
#'
#' @return a vector of variables corresponding to the projection.
#' @keywords internal

sample_vars_proj <- function(ids_x_na,
                             X,
                             projection_function = NULL,
                             normal_proj = TRUE) {
  if (is.null(projection_function)) {
    if (length(ids_x_na) == 1) {
      vars_na <- ids_x_na
    } else {
      n_var_na <- sample(seq_along(ids_x_na), size = 1)
      vars_na <- sample(ids_x_na, size = n_var_na, replace = FALSE)
    }

    vars_na <- sort(vars_na)
    vars <- vars_na

    if (normal_proj) {
      avail <- setdiff(seq_len(ncol(X)), vars_na)

      if (length(avail) > 0) {
        if (ncol(X) == 2) {
          dim_proj <- 1
        } else {
          dim_proj <- sample(seq_len(length(avail)), size = 1)
        }

        extra_vars <- if (length(avail) == 1) {
          avail
        } else {
          sample(avail, size = dim_proj, replace = FALSE)
        }

        vars <- c(vars_na, extra_vars)
      }
    }
  } else {
    vars <- c(ids_x_na, projection_function(X))
  }

  sort(unique(vars))
}

#' Merge singleton missingness patterns
#'
#' @description
#' Merges missingness patterns that occur only once (singleton patterns) into a
#' single pattern. If the merged pattern already exists among the current
#' patterns, the corresponding groups of observations are combined. Otherwise,
#' a new pattern is created and appended.
#'
#' @param patterns A numeric matrix where each row represents a unique
#' missingness pattern.
#' @param groups A list of integer vectors. Each element contains the indices of
#' observations corresponding to a given pattern in \code{patterns}.
#' @param ind_singletons An integer vector indicating indices of patterns in
#' \code{patterns} that occur only once.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{patterns}{Updated matrix of unique missingness patterns.}
#'   \item{groups}{Updated list of observation indices grouped by pattern.}
#' }
#' @keywords internal

merge_singleton_patterns <- function(patterns, groups, ind_singletons) {
  obs_to_merge <- unlist(groups[ind_singletons])

  pat_singletons <- ifelse(
    colSums(patterns[ind_singletons, ], na.rm = TRUE) == 0,
    0,
    1
  )

  same_pattern <- which(
    rowSums(
      patterns == rep(pat_singletons, each = nrow(patterns))
    ) == ncol(patterns)
  )

  if (length(same_pattern) > 0) {
    groups[[same_pattern]] <- unique(c(groups[[same_pattern]], obs_to_merge))
    groups <- groups[-setdiff(ind_singletons, same_pattern)]
    patterns <- patterns[-setdiff(ind_singletons, same_pattern), ]
  } else {
    groups <- groups[-ind_singletons]
    patterns <- patterns[-ind_singletons, ]

    groups <- c(groups, list(obs_to_merge))
    patterns <- rbind(patterns, pat_singletons)
  }

  list(patterns = patterns, groups = groups)
}

#' Extract and group missing-data patterns
#'
#' @description
#' Identifies unique missingness patterns in a data matrix and groups
#' observations according to these patterns. If more than one pattern occurs
#' only once, such singleton patterns are merged into a single group.
#'
#' @importFrom stats complete.cases
#'
#' @param X A matrix or data frame that may contain missing values.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{patterns}{A matrix of unique missingness patterns.}
#'   \item{groups}{A list of integer vectors giving row indices for each pattern.}
#'   \item{average_diff}{A logical indicating whether singleton patterns were merged.}
#' }
#'
#' @details
#' Missingness patterns are represented by a logical matrix obtained from
#' \code{is.na(X)}. Only rows containing at least one missing value are used
#' to define the unique patterns.
#'
#' If more than one pattern is represented by a single observation, these
#' singleton patterns are merged using \code{merge_singleton_patterns()}.
#'
#' @keywords internal

get_pattern_data <- function(X) {
  NA.pat <- is.na(X)
  NA_pat_unique <- unique(NA.pat[which(!complete.cases(X)), ])

  NA_pat_groups <- lapply(seq_len(nrow(NA_pat_unique)), function(i) {
    which(rowSums(abs(t(t(NA.pat) - NA_pat_unique[i, ]))) == 0)
  })

  lengths_groups <- lengths(NA_pat_groups)

  average_diff <- FALSE

  if (sum(lengths_groups == 1) > 1) {
    ind_singletons <- which(lengths_groups == 1)

    merged_res <- merge_singleton_patterns(
      patterns = NA_pat_unique,
      groups = NA_pat_groups,
      ind_singletons = ind_singletons
    )

    NA_pat_unique <- merged_res[["patterns"]]

    if (nrow(NA_pat_unique) == length(NA_pat_groups) - length(ind_singletons) + 1) {
      average_diff <- TRUE
    }

    NA_pat_groups <- merged_res[["groups"]]
  }

  list(
    patterns = NA_pat_unique,
    groups = NA_pat_groups,
    average_diff = average_diff
  )
}

#' Compute the imputation KL-based scoring rules
#'
#' @importFrom pbmcapply pbmclapply
#'
#' @inheritParams energy_Iscore_num
#'
#' @param imputation_func an imputing function. If \code{NULL}, please provide
#' imputed datasets \code{X_imp} and \code{m}.
#' @param X_imp a list of imputed datasets. If \code{NULL} it will be obtained
#' using \code{imputation_func}.
#' @param m the number of multiple imputations to consider, default to 5.
#' @param n_proj an integer specifying the number of projections to consider
#' for the score.
#' @param n_trees_per_proj an integer, the number of trees per projection.
#' @param min_node_size the minimum number of nodes in a tree.
#' @param n_cores an integer, the number of cores to use.
#' @param projection_function a function providing the user-specific projections.
#' @param ... used for compatibility
#'
#' @return  numeric value of the score obtained for provided imputation method.
#'
#' @examples
#' set.seed(111)
#' X <- random_mcar_data(100, 3, 0.2)
#' imputation_func <- exp_imputation
#' DR_IScore(X, imputation_func, m = 2, n_proj = 10, n_trees_per_proj = 2 )
#'
#'
#' @references
#' This method is described in detail in:
#'
#' Näf, Jeffrey, Meta-Lina Spohn, Loris Michel, and Nicolai Meinshausen. 2022.
#' “Imputation Scores.” https://arxiv.org/abs/2106.03742.
#'
#' @export
#'

DR_IScore <- function(X,
                      imputation_func = NULL,
                      X_imp = NULL,
                      m = 5,
                      n_proj = 100,
                      n_trees_per_proj = 5,
                      min_node_size = 10,
                      n_cores = 1,
                      projection_function = NULL,
                      ...) {
  if (is.null(X_imp) & is.null(imputation_func)) {
    stop("You must provide one of imputation_func or X_imp!")
  }

  if (is.null(X_imp)) {
    X_imp <- lapply(seq_len(m), function(i) {
      X_imp <- try({
        imputation_func(X)
      })

      if (inherits(X_imp, "try-error") | any(is.na(X_imp))) {
        stop("Errored imputing X using provided imputation_func!")
      }

      X_imp
    })
  }

  X <- as.matrix(X)

  pattern_data <- get_pattern_data(X)

  NA_pat_unique <- pattern_data[["patterns"]]
  NA_pat_groups <- pattern_data[["groups"]]
  average_diff <- pattern_data[["average_diff"]]

  res <- pbmcapply::pbmclapply(seq_len(nrow(NA_pat_unique)), function(j) {
    is_last_and_merged <- j == nrow(NA_pat_unique) & average_diff

    group_j <- NA_pat_groups[[j]]
    all_ids <- seq_len(nrow(X))

    is_singleton <- length(group_j) == 1

    normal_proj <- ifelse(is_last_and_merged || is_singleton, FALSE, TRUE)
    parts <- c(1, 2)

    if (is_last_and_merged | is_singleton) {
      parts <- 1
    }

    scores.all <- lapply(parts, function(part) {
      if (is_last_and_merged || is_singleton) {
        ids_pat_test <- group_j
        ids_pat_train <- all_ids
      } else {
        half <- floor(length(group_j) / 2)

        if (part == 1) {
          ids_pat_test <- group_j[seq_len(half)]
        } else {
          ids_pat_test <- group_j[-seq_len(half)]
        }

        ids_pat_train <- setdiff(all_ids, ids_pat_test)
      }

      scores <- sapply(seq_along(X_imp), function(set) {
        object_dr <- densityRatioScore(
          X = X[ids_pat_train, ],
          X_imp = X_imp[[set]][ids_pat_train, ],
          pattern = NA_pat_unique[j, ],
          n_proj = n_proj,
          n_trees_per_proj = n_trees_per_proj,
          min_node_size = min_node_size,
          projection_function = projection_function,
          normal_proj = normal_proj
        )

        if (any(is.na(object_dr)) || any(is.null(object_dr))) {
          score_dr_kl <- NA
        } else {
          Z <- unname(as.matrix(X_imp[[set]])[ids_pat_test, , drop = FALSE])

          score_dr_kl <- compute_drScore(
            object = object_dr,
            Z = Z,
            n_trees_per_proj = n_trees_per_proj,
            n_proj = n_proj
          )

          if (is_last_and_merged) {
            score_dr_kl <- unlist(score_dr_kl)
          } else {
            score_dr_kl <- mean(score_dr_kl)
          }
        }

        score_dr_kl
      })

      scores
    })

    mean(unlist(scores.all), na.rm = TRUE)

  }, mc.cores = n_cores)

  mean(unlist(res), na.rm = TRUE)
}
