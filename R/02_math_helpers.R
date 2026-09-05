# Basic mathematical helpers
#
# This file follows CONVENTIONS.md. Names from the pasted code have been
# converted from period separators or unexplained capitalization to snake_case.
# Compiled accelerators are registered when the package is installed.

# Return the most frequent value in x, using the first value to break ties.
#
# If na_rm is FALSE, NA is treated as a possible modal value. If removal leaves
# no observations, return a typed NA value.
#' Compute the statistical mode
#'
#' Returns the most frequent value in `x`, using the first value to break ties.
#' If `na_rm = FALSE`, `NA` is a possible modal value; if removal leaves no
#' observations, a typed missing value is returned.
#' @param x Input vector.
#' @param na_rm Remove missing values before finding the mode.
#' @return A scalar with the same basic type as `x`.
#' @examples modal(c(1, 2, 2, 3))
#' @seealso [nmi()]
#' @name modal
#' @export
modal <- function(x, na_rm = FALSE) {
  if (length(na_rm) != 1L || !is.logical(na_rm) || is.na(na_rm)) {
    stop("na_rm must be TRUE or FALSE.", call. = FALSE)
  }
  if (na_rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) == 0L) {
    return(x[NA_integer_][1L])
  }

  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

# Validate inputs shared by vector error/loss functions.
validate_error_inputs <- function(x, y, na_rm) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric.", call. = FALSE)
  }
  if (length(x) != length(y)) {
    stop("x and y must have equal lengths.", call. = FALSE)
  }
  if (length(na_rm) != 1L || !is.logical(na_rm) || is.na(na_rm)) {
    stop("na_rm must be TRUE or FALSE.", call. = FALSE)
  }

  invisible(TRUE)
}

# Compute the sum of squared errors between equal-length numeric vectors.
#' Compute squared-error losses
#'
#' `sse()` returns the sum of squared errors; `mse()` returns their mean.
#' @param x,y Equal-length numeric vectors.
#' @param na_rm Remove missing comparisons.
#' @param validate_inputs Validate input types, lengths, and `na_rm`.
#' @return A nonnegative numeric scalar.
#' @examples
#' sse(c(1, 2), c(1, 3))
#' mse(c(1, 2), c(1, 3))
#' @seealso [sae()], [mae()]
#' @name squared_error_losses
#' @export
sse <- function(x, y, na_rm = FALSE, validate_inputs = TRUE) {
  if (validate_inputs) validate_error_inputs(x, y, na_rm)

  sum((x - y)^2, na.rm = na_rm)
}

# Compute the sum of absolute errors between equal-length numeric vectors.
#' Compute absolute-error losses
#'
#' `sae()` returns the sum of absolute errors; `mae()` returns their mean.
#' @param x,y Equal-length numeric vectors.
#' @param na_rm Remove missing comparisons.
#' @param validate_inputs Validate input types, lengths, and `na_rm`.
#' @return A nonnegative numeric scalar.
#' @examples
#' sae(c(1, 2), c(1, 3))
#' mae(c(1, 2), c(1, 3))
#' @seealso [sse()], [mse()]
#' @name absolute_error_losses
#' @export
sae <- function(x, y, na_rm = FALSE, validate_inputs = TRUE) {
  if (validate_inputs) validate_error_inputs(x, y, na_rm)

  sum(abs(x - y), na.rm = na_rm)
}

# Compute the mean squared error between equal-length numeric vectors.
#' @rdname squared_error_losses
#' @export
mse <- function(x, y, na_rm = FALSE, validate_inputs = TRUE) {
  if (validate_inputs) validate_error_inputs(x, y, na_rm)

  mean((x - y)^2, na.rm = na_rm)
}

# Compute the mean absolute error between equal-length numeric vectors.
#' @rdname absolute_error_losses
#' @export
mae <- function(x, y, na_rm = FALSE, validate_inputs = TRUE) {
  if (validate_inputs) validate_error_inputs(x, y, na_rm)

  mean(abs(x - y), na.rm = na_rm)
}

# Compute empirical mutual information normalized by joint entropy.
#
# This preserves the requested normalization MI(g_1, g_2) / H(g_1, g_2),
# which differs from normalizations based on the two marginal entropies. Use
# entropy when it is installed; otherwise compute the same empirical plug-in
# quantities directly. A constant joint labeling has zero entropy, making the
# normalized value undefined, and returns NaN.
#' Compute normalized mutual information
#'
#' Computes empirical mutual information divided by joint entropy. This
#' normalization differs from definitions based on marginal entropies. A
#' constant joint labeling has zero entropy and returns `NaN`.
#' @param g_1,g_2 Equal-length, non-missing atomic label vectors.
#' @return A numeric agreement score.
#' @examples nmi(c(1, 1, 2, 2), c(2, 2, 1, 1))
#' @seealso [label_match_greedy()], [pair_nmi_loss()]
#' @name nmi
#' @export
nmi <- function(g_1, g_2) {
  if (!is.atomic(g_1) || !is.null(dim(g_1)) ||
      !is.atomic(g_2) || !is.null(dim(g_2))) {
    stop("g_1 and g_2 must be atomic vectors.", call. = FALSE)
  }
  if (length(g_1) != length(g_2)) {
    stop("g_1 and g_2 must have equal lengths.", call. = FALSE)
  }
  if (length(g_1) == 0L) {
    stop("g_1 and g_2 must not be empty.", call. = FALSE)
  }
  if (anyNA(g_1) || anyNA(g_2)) {
    stop("g_1 and g_2 must not contain missing values.", call. = FALSE)
  }

  cross_table <- table(g_1, g_2)
  if (requireNamespace("entropy", quietly = TRUE)) {
    mutual_information <- entropy::mi.empirical(
      as.matrix(cross_table)
    )
    joint_entropy <- entropy::entropy.empirical(
      as.vector(cross_table)
    )
  } else {
    joint_probabilities <- as.numeric(cross_table) / sum(cross_table)
    positive_joint <- joint_probabilities > 0
    joint_entropy <- -sum(
      joint_probabilities[positive_joint] *
        log(joint_probabilities[positive_joint])
    )
    row_probabilities <- rowSums(cross_table) / sum(cross_table)
    column_probabilities <- colSums(cross_table) / sum(cross_table)
    expected_probabilities <- tcrossprod(
      row_probabilities,
      column_probabilities
    )
    mutual_information <- sum(
      joint_probabilities[positive_joint] * log(
        joint_probabilities[positive_joint] /
          expected_probabilities[positive_joint]
      )
    )
  }

  if (joint_entropy == 0) return(NaN)
  as.numeric(mutual_information / joint_entropy)
}

# Clip every value in a numeric vector or matrix to inclusive lower/upper bounds.
#
# Missing values remain missing. Infinite values are replaced when the
# corresponding bound is finite. Names, dimensions, and dimension names are
# preserved by the elementwise pmax/pmin operations.
#' Clip numeric values to an interval
#'
#' Restricts each value to inclusive bounds while preserving missing values,
#' names, dimensions, and dimension names.
#' @param x Numeric vector or matrix.
#' @param lower_clip,upper_clip Inclusive clipping bounds.
#' @return An object shaped like `x` with clipped values.
#' @examples clip_values(c(-1, 0.5, 2), 0, 1)
#' @seealso [clip_probabilities()]
#' @name clip_values
#' @export
clip_values <- function(x, lower_clip = -Inf, upper_clip = Inf) {
  if (!is.numeric(x)) {
    stop("x must be a numeric vector or matrix.", call. = FALSE)
  }
  if (length(dim(x)) > 2L) {
    stop("x must be a vector or matrix, not a higher-dimensional array.",
         call. = FALSE)
  }
  if (length(lower_clip) != 1L ||
      !is.numeric(lower_clip) ||
      is.na(lower_clip)) {
    stop("lower_clip must be one non-missing numeric value.", call. = FALSE)
  }
  if (length(upper_clip) != 1L ||
      !is.numeric(upper_clip) ||
      is.na(upper_clip)) {
    stop("upper_clip must be one non-missing numeric value.", call. = FALSE)
  }
  if (lower_clip > upper_clip) {
    stop("lower_clip must not exceed upper_clip.", call. = FALSE)
  }

  pmin(pmax(x, lower_clip), upper_clip)
}

# Clip probabilities to [eps, 1 - eps] to protect logs and inverse links.
#' Clip probability values safely
#'
#' Restricts probabilities to `[eps, 1 - eps]` to protect logarithms and
#' inverse-link calculations.
#' @param P Numeric probability vector or matrix.
#' @param eps Finite clipping constant strictly between zero and one half.
#' @return An object shaped like `P` with clipped probabilities.
#' @examples clip_probabilities(c(0, 0.5, 1), eps = 0.01)
#' @seealso [clip_values()], [bin_dev()]
#' @name clip_probabilities
#' @export
clip_probabilities <- function(P, eps = 1e-6) {
  if (length(eps) != 1L ||
      !is.numeric(eps) ||
      is.na(eps) ||
      !is.finite(eps) ||
      eps <= 0 ||
      eps >= 0.5) {
    stop("eps must be one finite number strictly between 0 and 0.5.",
         call. = FALSE)
  }

  clip_values(P, lower_clip = eps, upper_clip = 1 - eps)
}

# Compute binary cross-entropy loss for binary responses and probabilities.
#
# Probabilities are clipped to [epsilon, 1 - epsilon] before taking logarithms.
# Despite the historical name, this returns the summed negative log-likelihood,
# without the factor of two sometimes included in definitions of deviance.
#' Compute binary deviance losses
#'
#' `bin_dev()` returns summed binary cross-entropy after clipping predictions;
#' despite its historical name it does not multiply by two.
#' `bin_dev_mean()` averages over retained response-probability pairs.
#' @param x Binary response vector.
#' @param y Numeric predicted-probability vector.
#' @param epsilon Finite clipping constant strictly between zero and one half.
#' @param na_rm Remove comparisons containing missing values.
#' @param validate_inputs Validate types, lengths, responses, and predictions.
#' @return A nonnegative sum from `bin_dev()` or mean from `bin_dev_mean()`.
#' @examples
#' bin_dev(c(0, 1), c(0.1, 0.8))
#' bin_dev_mean(c(0, 1), c(0.1, 0.8))
#' @seealso [auc()], [clip_probabilities()]
#' @name binary_deviance_losses
#' @export
bin_dev <- function(
    x, y, epsilon = 1e-5, na_rm = TRUE, validate_inputs = TRUE) {
  if (validate_inputs && (!is.numeric(x) || !is.numeric(y))) {
    stop("x and y must be numeric.", call. = FALSE)
  }
  if (validate_inputs && length(x) != length(y)) {
    stop("x and y must have equal lengths.", call. = FALSE)
  }
  if (length(epsilon) != 1L ||
      !is.numeric(epsilon) ||
      is.na(epsilon) ||
      !is.finite(epsilon) ||
      epsilon <= 0 ||
      epsilon >= 0.5) {
    stop("epsilon must be one finite number strictly between 0 and 0.5.",
         call. = FALSE)
  }
  if (length(na_rm) != 1L || !is.logical(na_rm) || is.na(na_rm)) {
    stop("na_rm must be TRUE or FALSE.", call. = FALSE)
  }

  if (validate_inputs) {
    observed_x <- x[!is.na(x)]
    if (any(!is.finite(observed_x)) || any(!observed_x %in% c(0, 1))) {
      stop("Non-missing values of x must be binary (0 or 1).", call. = FALSE)
    }
    observed_y <- y[!is.na(y)]
    if (any(!is.finite(observed_y))) {
      stop("Non-missing values of y must be finite.", call. = FALSE)
    }
  }

  y_clipped <- clip_probabilities(y, eps = epsilon)
  losses <- -x * log(y_clipped) - (1 - x) * log1p(-y_clipped)
  sum(losses, na.rm = na_rm)
}

# Compute mean binary cross-entropy over the observed response/probability pairs.
#' @rdname binary_deviance_losses
#' @export
bin_dev_mean <- function(
    x, y, epsilon = 1e-5, na_rm = TRUE, validate_inputs = TRUE) {
  total_loss <- bin_dev(
    x, y, epsilon = epsilon, na_rm = na_rm,
    validate_inputs = validate_inputs
  )
  if (!na_rm) {
    return(total_loss / length(x))
  }

  observed <- !is.na(x) & !is.na(y)
  total_loss / sum(observed)
}

# Compute binary ROC AUC from labels and prediction scores.
#
# Use the compiled auc_cpp implementation when it is loaded and working;
# otherwise, use an equivalent base-R rank calculation. This makes the helper
# usable when Rcpp, a compiler, the shared library, or the exported function is
# unavailable. Ties receive their average rank in both implementations.
#' Compute AUC scores and losses
#'
#' `auc()` computes binary ROC AUC using average ranks for tied scores and a
#' compiled implementation when available. `auc_as_loss()` returns one minus
#' AUC so smaller values are better.
#' @param A Binary response vector.
#' @param P Numeric prediction-score vector.
#' @param use_cpp Try the registered compiled implementation before the R
#'   implementation.
#' @param validate_inputs Validate binary responses and finite scores.
#' @return A numeric scalar containing AUC or one minus AUC.
#' @examples
#' auc(c(0, 1, 1), c(0.1, 0.7, 0.8), use_cpp = FALSE)
#' auc_as_loss(c(0, 1, 1), c(0.1, 0.7, 0.8), use_cpp = FALSE)
#' @seealso [bin_dev()], [bin_dev_mean()]
#' @name auc_losses
#' @export
auc <- function(A, P, use_cpp = TRUE, validate_inputs = TRUE) {
  group <- as.numeric(A)
  predictions <- as.numeric(P)
  if (validate_inputs && (!is.numeric(group) || !is.numeric(predictions))) {
    stop("group and predictions must be numeric.", call. = FALSE)
  }
  if (validate_inputs && length(group) != length(predictions)) {
    stop("group and predictions must have equal lengths.", call. = FALSE)
  }
  if (validate_inputs && any(!is.finite(predictions))) {
    stop("predictions must be finite.", call. = FALSE)
  }
  if (validate_inputs &&
      (any(!is.finite(group)) || any(!group %in% c(0, 1)))) {
    stop("group must contain only 0 and 1.", call. = FALSE)
  }
  if (length(use_cpp) != 1L || !is.logical(use_cpp) || is.na(use_cpp)) {
    stop("use_cpp must be TRUE or FALSE.", call. = FALSE)
  }

  if (use_cpp && exists("auc_cpp", mode = "function", inherits = TRUE, envir = environment())) {
    compiled_auc <- get("auc_cpp", mode = "function", inherits = TRUE, envir = environment())
    compiled_result <- tryCatch(
      list(
        success = TRUE,
        value = compiled_auc(
          group = as.numeric(group),
          predictions = as.numeric(predictions)
        )
      ),
      error = function(e) {
        message(
          "auc_cpp failed; using the base-R AUC implementation: ",
          conditionMessage(e)
        )
        list(success = FALSE, value = NULL)
      }
    )
    if (compiled_result$success) {
      return(compiled_result$value)
    }
  }

  n1 <- sum(group == 1)
  n0 <- length(group) - n1
  if (n0 == 0L || n1 == 0L) {
    return(NaN)
  }

  prediction_ranks <- rank(predictions, ties.method = "average")
  positive_rank_sum <- sum(prediction_ranks[group == 1])
  u_statistic <- positive_rank_sum - n1 * (n1 + 1) / 2
  as.numeric(u_statistic / (as.double(n0) * n1))
}

# Convert ROC AUC, where larger is better, into a loss where smaller is better.
#' @rdname auc_losses
#' @export
auc_as_loss <- function(A, P, use_cpp = TRUE, validate_inputs = TRUE) {
  1 - auc(
    A = A, P = P, use_cpp = use_cpp,
    validate_inputs = validate_inputs
  )
}

# Compute log(1 + exp(x)) without overflow or avoidable underflow.
#' Compute the softplus transform
#'
#' Evaluates the numerically stable softplus transform elementwise.
#' @param x Numeric vector or matrix.
#' @return A numeric object matching the shape of `x`.
#' @examples softplus(c(-1, 0, 1))
#' @seealso [sigmoid()]
#' @name softplus
#' @export
softplus <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# Compute the logistic sigmoid with a numerically stable base-R implementation.
#' Compute the sigmoid transform
#'
#' Evaluates the numerically stable logistic sigmoid elementwise.
#' @param x Numeric vector or matrix.
#' @return A numeric object matching the shape of `x`.
#' @examples sigmoid(c(-1, 0, 1))
#' @seealso [softplus()]
#' @name sigmoid
#' @export
sigmoid <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  stats::plogis(x)
}

# Validate inputs shared by soft/hard extrema-selection functions.
validate_extrema_inputs <- function(x, na_rm, temperature = NULL) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("x must be a numeric vector.", call. = FALSE)
  }
  if (length(na_rm) != 1L || !is.logical(na_rm) || is.na(na_rm)) {
    stop("na_rm must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(temperature) &&
      (length(temperature) != 1L ||
       !is.numeric(temperature) ||
       is.na(temperature) ||
       !is.finite(temperature) ||
       temperature <= 0)) {
    stop("temperature must be one finite positive number.", call. = FALSE)
  }

  invisible(TRUE)
}

# Convert scores to normalized exponential weights favoring larger values.
#
# Missing positions receive zero weight when na_rm is TRUE. If na_rm is FALSE,
# any missing score produces an all-NA result. Positive infinities split all
# mass equally; an all-negative-infinity vector is treated as a complete tie.
#' Compute hard and soft extrema
#'
#' `softmax()` and `softmin()` return normalized exponential weights.
#' `hardmax()` and `hardmin()` return zero-one selectors marking all ties.
#' The `which_soft*()` functions sample one index from the soft weights; the
#' `which_hard*()` functions return every tied hard-extremum index.
#' @param x Numeric vector whose extrema are selected.
#' @param temperature Positive soft-selection temperature.
#' @param na_rm Remove missing values before selection.
#' @return A numeric weight vector or selected index vector, according to the
#'   function called.
#' @examples
#' softmax(c(1, 2, 3))
#' hardmin(c(2, 1, 1))
#' which_hardmax(c(1, 3, 3))
#' @seealso [softplus()], [sigmoid()]
#' @name extrema_selectors
#' @export
softmax <- function(x, temperature = 1, na_rm = FALSE) {
  validate_extrema_inputs(x, na_rm, temperature)

  weights <- rep(0, length(x))
  names(weights) <- names(x)
  observed <- !is.na(x)
  if (!na_rm && any(!observed)) {
    weights[] <- NA_real_
    return(weights)
  }
  if (!any(observed)) {
    weights[] <- NA_real_
    return(weights)
  }

  observed_x <- x[observed]
  positive_infinite <- is.infinite(observed_x) & observed_x > 0
  if (any(positive_infinite)) {
    observed_weights <- as.numeric(positive_infinite) / sum(positive_infinite)
  } else if (all(is.infinite(observed_x) & observed_x < 0)) {
    observed_weights <- rep(1 / length(observed_x), length(observed_x))
  } else {
    shifted_x <- (observed_x - max(observed_x)) / temperature
    observed_weights <- exp(shifted_x)
    observed_weights <- observed_weights / sum(observed_weights)
  }

  weights[observed] <- observed_weights
  weights
}

# Convert scores to normalized exponential weights favoring smaller values.
#' @rdname extrema_selectors
#' @export
softmin <- function(x, temperature = 1, na_rm = FALSE) {
  validate_extrema_inputs(x, na_rm, temperature)
  softmax(-x, temperature = temperature, na_rm = na_rm)
}

# Return a 0-1 selector marking every maximal value in x.
#
# Tied maxima all receive one. Missing-value behavior matches softmax.
#' @rdname extrema_selectors
#' @export
hardmax <- function(x, na_rm = FALSE) {
  validate_extrema_inputs(x, na_rm)

  weights <- rep(0, length(x))
  names(weights) <- names(x)
  observed <- !is.na(x)
  if (!na_rm && any(!observed)) {
    weights[] <- NA_real_
    return(weights)
  }
  if (!any(observed)) {
    weights[] <- NA_real_
    return(weights)
  }

  maximum <- max(x[observed])
  weights[observed & x == maximum] <- 1
  weights
}

# Return a 0-1 selector marking every minimal value in x.
#' @rdname extrema_selectors
#' @export
hardmin <- function(x, na_rm = FALSE) {
  validate_extrema_inputs(x, na_rm)
  hardmax(-x, na_rm = na_rm)
}

# Sample one index according to the softmax weights.
#' @rdname extrema_selectors
#' @export
which_softmax <- function(x, temperature = 1, na_rm = FALSE) {
  probabilities <- softmax(x, temperature = temperature, na_rm = na_rm)
  if (length(probabilities) == 0L || all(is.na(probabilities))) {
    return(NA_integer_)
  }
  sample.int(length(probabilities), size = 1L, prob = probabilities)
}

# Sample one index according to the softmin weights.
#' @rdname extrema_selectors
#' @export
which_softmin <- function(x, temperature = 1, na_rm = FALSE) {
  probabilities <- softmin(x, temperature = temperature, na_rm = na_rm)
  if (length(probabilities) == 0L || all(is.na(probabilities))) {
    return(NA_integer_)
  }
  sample.int(length(probabilities), size = 1L, prob = probabilities)
}

# Return all indices tied for the maximum value.
#' @rdname extrema_selectors
#' @export
which_hardmax <- function(x, na_rm = FALSE) {
  weights <- hardmax(x, na_rm = na_rm)
  if (length(weights) == 0L || all(is.na(weights))) {
    return(NA_integer_)
  }
  which(weights == 1)
}

# Return all indices tied for the minimum value.
#' @rdname extrema_selectors
#' @export
which_hardmin <- function(x, na_rm = FALSE) {
  weights <- hardmin(x, na_rm = na_rm)
  if (length(weights) == 0L || all(is.na(weights))) {
    return(NA_integer_)
  }
  which(weights == 1)
}

# Form the additive outer combination with y as rows and x as columns.
#
# The compiled implementation is optional. If it is unavailable or fails, the
# wrapper uses base::outer with identical orientation and output semantics.
#' Add vectors by outer expansion
#'
#' Forms all pairwise sums, using a registered compiled implementation when
#' requested and available and otherwise [base::outer()].
#' @param x,y Numeric vectors.
#' @param operator Outer operation; currently only `"+"` is supported.
#' @param use_cpp Try the compiled implementation before the R implementation.
#' @return A numeric matrix with `length(x)` rows and `length(y)` columns.
#' @examples outer_add(1:2, 3:4, use_cpp = FALSE)
#' @seealso [procrustes()]
#' @name outer_add
#' @export
outer_add <- function(x, y, operator = "+", use_cpp = TRUE) {
  if (!identical(operator, "+")) {
    stop("outer_add currently supports only operator = '+'.", call. = FALSE)
  }
  if (!is.numeric(x) || !is.null(dim(x)) ||
      !is.numeric(y) || !is.null(dim(y))) {
    stop("x and y must be numeric vectors.", call. = FALSE)
  }
  if (length(use_cpp) != 1L || !is.logical(use_cpp) || is.na(use_cpp)) {
    stop("use_cpp must be TRUE or FALSE.", call. = FALSE)
  }

  result <- NULL
  if (use_cpp && exists("outer_add_cpp", mode = "function", inherits = TRUE, envir = environment())) {
    compiled_outer_add <- get(
      "outer_add_cpp",
      mode = "function",
      inherits = TRUE, envir = environment()
    )
    compiled_result <- tryCatch(
      list(
        success = TRUE,
        value = compiled_outer_add(
          x = as.numeric(x),
          y = as.numeric(y)
        )
      ),
      error = function(e) {
        message(
          "outer_add_cpp failed; using base::outer: ",
          conditionMessage(e)
        )
        list(success = FALSE, value = NULL)
      }
    )
    if (compiled_result$success) {
      result <- compiled_result$value
    }
  }
  if (is.null(result)) {
    result <- base::outer(y, x, FUN = "+")
  }

  dimnames(result) <- list(names(y), names(x))
  result
}

# Align X to X_star with an orthogonal Procrustes transformation.
#
# Translation-only alignment uses the optional compiled implementation when it
# is available. Every other path uses BLAS-backed R matrix operations. X may be
# padded with zero columns when X_star has a larger embedding dimension.
#' Align point configurations by Procrustes transformation
#'
#' Aligns `X` to `X_star` by an orthogonal transformation with optional
#' centering and dilation. Inputs may have different column counts; the
#' narrower configuration is padded with zero columns. A compiled translated,
#' non-dilated path is used when requested and available.
#' @param X Numeric point-configuration matrix to align.
#' @param X_star Numeric reference point-configuration matrix with the same
#'   number of rows as `X`.
#' @param translate Center configurations and estimate translation.
#' @param dilate Estimate and apply a dilation factor.
#' @param sumsq Include residual sum of squares.
#' @param use_cpp Try the registered compiled implementation.
#' @param validate_inputs Validate matrices and logical controls.
#' @return A list containing aligned `X_new` and `rotation`, with optional
#'   `translation`, `scale`, and `sse` components.
#' @examples
#' X <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
#' procrustes(X, X, use_cpp = FALSE)$X_new
#' @seealso [ase()], [generate_latent_positions()]
#' @name procrustes
#' @export
procrustes <- function(
    X,
    X_star,
    translate = FALSE,
    dilate = FALSE,
    sumsq = FALSE,
    use_cpp = TRUE,
    validate_inputs = TRUE) {
  logical_inputs <- list(
    translate = translate,
    dilate = dilate,
    sumsq = sumsq,
    use_cpp = use_cpp,
    validate_inputs = validate_inputs
  )
  invalid_logical <- vapply(
    logical_inputs,
    function(value) {
      length(value) != 1L || !is.logical(value) || is.na(value)
    },
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_inputs)[invalid_logical], collapse = ", "),
      " must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  if (validate_inputs) {
    if (is.null(dim(X)) || length(dim(X)) != 2L ||
        is.null(dim(X_star)) || length(dim(X_star)) != 2L) {
      stop("X and X_star must be matrix-like objects.", call. = FALSE)
    }
    if (!(is.numeric(X) || inherits(X, "Matrix")) ||
        !(is.numeric(X_star) || inherits(X_star, "Matrix"))) {
      stop("X and X_star must be numeric.", call. = FALSE)
    }
  }

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(X_star)) X_star <- as.matrix(X_star)
  if (validate_inputs &&
      (any(!is.finite(X)) || any(!is.finite(X_star)))) {
    stop("X and X_star must contain only finite values.", call. = FALSE)
  }
  N <- nrow(X)
  if (validate_inputs && (N < 1L || N != nrow(X_star))) {
    stop("X and X_star must have the same positive number of rows.",
         call. = FALSE)
  }
  P <- ncol(X)
  P_star <- ncol(X_star)
  if (validate_inputs && P_star < 1L) {
    stop("X_star must contain at least one column.", call. = FALSE)
  }
  if (validate_inputs && P > P_star) {
    stop("X cannot have more columns than X_star.", call. = FALSE)
  }
  if (P < P_star) {
    warning(
      "X was padded with zeros to match the columns of X_star.",
      call. = FALSE,
      immediate. = TRUE
    )
    X <- cbind(X, matrix(0, nrow = N, ncol = P_star - P))
  }

  native_result <- NULL
  if (translate && !dilate && use_cpp &&
      exists("procrustes_translated_cpp", mode = "function", inherits = TRUE, envir = environment())) {
    compiled_procrustes <- get(
      "procrustes_translated_cpp",
      mode = "function",
      inherits = TRUE, envir = environment()
    )
    compiled_result <- tryCatch(
      list(
        success = TRUE,
        value = compiled_procrustes(X = X, X_star = X_star)
      ),
      error = function(e) {
        message(
          "procrustes_translated_cpp failed; using the R implementation: ",
          conditionMessage(e)
        )
        list(success = FALSE, value = NULL)
      }
    )
    if (compiled_result$success) {
      native_result <- compiled_result$value
    }
  }

  if (!is.null(native_result)) {
    X_new <- unname(native_result$X_new)
    rotation <- native_result$rotation
    translation <- as.numeric(native_result$translation)
    scale <- 1
  } else {
    if (translate) {
      mean_x <- colMeans(X)
      mean_x_star <- colMeans(X_star)
      cross <- crossprod(X_star, X) -
        N * tcrossprod(mean_x_star, mean_x)
    } else {
      mean_x <- NULL
      mean_x_star <- NULL
      cross <- crossprod(X_star, X)
    }

    decomposition <- svd(cross)
    rotation <- tcrossprod(decomposition$v, decomposition$u)

    if (dilate) {
      denominator <- if (translate) {
        sum(X^2) - N * sum(mean_x^2)
      } else {
        sum(X^2)
      }
      denominator_tolerance <- .Machine$double.eps * max(1, sum(X^2))
      if (!is.finite(denominator) || denominator <= denominator_tolerance) {
        stop("Cannot dilate because centered X has zero numerical variation.",
             call. = FALSE)
      }
      scale <- sum(cross * t(rotation)) / denominator
    } else {
      scale <- 1
    }

    scaled_rotated <- scale * tcrossprod(X, t(rotation))
    translation <- if (translate) {
      mean_x_star - colMeans(scaled_rotated)
    } else {
      0
    }
    X_new <- unname(if (translate) {
      sweep(scaled_rotated, 2L, translation, `+`)
    } else {
      scaled_rotated
    })
  }

  dimnames(X_new) <- dimnames(X_star)
  result <- list(X_new = X_new, rotation = rotation)
  if (translate) {
    result$translation <- translation
  }
  if (dilate) {
    result$scale <- scale
  }
  if (sumsq) {
    result$sse <- sum((X_star - X_new)^2)
  }
  result
}

# Return the fraction of entries in a matrix-like object that are nonzero.
#
# This is the edge/entry density called "sparsity" by the decomposition APIs.
# The explicit name avoids ambiguity with the fraction of zero entries.
#' Compute matrix density
#'
#' Returns the fraction of matrix entries that are nonzero. This is the entry
#' density called `"sparsity"` by the decomposition APIs.
#' @param A A finite, nonempty dense or sparse matrix.
#' @return A numeric scalar between zero and one.
#' @examples matrix_density(diag(3))
#' @seealso [generate_adjacency()]
#' @name matrix_density
#' @export
matrix_density <- function(A) {
  if (is.null(dim(A)) || length(dim(A)) != 2L) {
    stop("A must be a matrix-like object.", call. = FALSE)
  }
  if (any(dim(A) < 1L)) {
    stop("A must have positive dimensions.", call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }

  as.numeric(sum(A != 0)) / (as.double(nrow(A)) * ncol(A))
}

# Test matrix symmetry for both base and Matrix-package matrix classes.
#' Test matrix symmetry internally
#'
#' Tests base and Matrix-package matrix classes for numerical symmetry.
#' @param A A dense or sparse matrix.
#' @return A logical scalar.
#' @keywords internal
#' @name is_symmetric_matrix
is_symmetric_matrix <- function(A) {
  if (inherits(A, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The Matrix package is required for sparse matrix inputs.",
           call. = FALSE)
    }
    return(isTRUE(Matrix::isSymmetric(A)))
  }
  # Node names need not be duplicated on both dimensions for numerical
  # symmetry; decomposition routines care about values, not dimname layout.
  isTRUE(isSymmetric(A, check.attributes = FALSE))
}

# Construct an unnormalized or symmetric normalized graph Laplacian.
#
# When tau is positive, regularize the adjacency before forming the Laplacian:
# A_tau = A + tau * mean(deg) / n. Thus tau is dimensionless, tau = 1 adds
# one average degree across every row, and tau = 0 recovers the usual matrix.
#' Construct a graph Laplacian
#'
#' Constructs an unnormalized or symmetric normalized graph Laplacian. For
#' positive `tau`, the adjacency is first regularized by adding
#' `tau * mean(degree) / n` to every entry; `tau = 0` gives the usual matrix.
#' @param A A finite symmetric square dense or sparse adjacency matrix.
#' @param normalized Construct the symmetric normalized Laplacian.
#' @param tau Nonnegative dimensionless regularization strength.
#' @return A dense or sparse Laplacian matrix.
#' @examples
#' A <- matrix(c(0, 1, 1, 0), 2, 2)
#' graph_laplacian(A, normalized = FALSE)
#' @seealso [spectral_cluster()], [ase()]
#' @name graph_laplacian
#' @export
graph_laplacian <- function(A, normalized = TRUE, tau = 0) {
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
    stop("A must be a square matrix-like object.", call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric.", call. = FALSE)
  }
  if (length(normalized) != 1L ||
      !is.logical(normalized) ||
      is.na(normalized)) {
    stop("normalized must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(tau) != 1L ||
      !is.numeric(tau) ||
      is.na(tau) ||
      !is.finite(tau) ||
      tau < 0) {
    stop("tau must be one finite non-negative number.", call. = FALSE)
  }

  is_sparse <- inherits(A, "Matrix")
  if (is_sparse && !requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for sparse matrix inputs.",
         call. = FALSE)
  }
  deg <- if (is_sparse) Matrix::rowSums(A) else rowSums(A)
  if (tau > 0) {
    regularization <- tau * mean(deg) / nrow(A)
    A <- A + regularization
    # Positive all-pairs regularization is dense even when A is sparse.
    # Matrix::rowSums and Matrix::diag are also required for the dense
    # dgeMatrix produced by adding a scalar to a sparse Matrix object.
    is_sparse <- inherits(A, "Matrix")
    deg <- if (is_sparse) Matrix::rowSums(A) else rowSums(A)
  }
  if (normalized && any(deg < 0)) {
    stop("Normalized Laplacians require non-negative node degrees.",
         call. = FALSE)
  }

  if (!normalized) {
    L <- -A
    if (is_sparse) {
      Matrix::diag(L) <- Matrix::diag(L) + as.numeric(deg)
    } else {
      diag(L) <- diag(L) + as.numeric(deg)
    }
    return(L)
  }

  inverse_sqrt_deg <- numeric(length(deg))
  positive_deg <- deg > 0
  inverse_sqrt_deg[positive_deg] <- 1 / sqrt(deg[positive_deg])
  normalized_A <- if (is_sparse) {
    scaling_matrix <- Matrix::Diagonal(x = inverse_sqrt_deg)
    crossprod(scaling_matrix, tcrossprod(A, scaling_matrix))
  } else {
    row_scaled_A <- sweep(A, 1L, inverse_sqrt_deg, `*`)
    sweep(row_scaled_A, 2L, inverse_sqrt_deg, `*`)
  }
  L <- -normalized_A
  if (is_sparse) {
    Matrix::diag(L) <- Matrix::diag(L) + as.numeric(positive_deg)
  } else {
    diag(L) <- diag(L) + as.numeric(positive_deg)
  }
  L
}

# Compute selected eigenpairs of a finite symmetric matrix.
#
# The partial RSpectra engine falls back to base::eigen when unavailable or
# unsuccessful unless force_engine is TRUE. The irlba option is retained for
# input compatibility but is redirected because an SVD is not a signed
# eigendecomposition. `scale_by = "sparsity"` uses matrix_density(A).
#' Compute an eigendecomposition
#'
#' Computes selected eigenpairs of a finite symmetric matrix. Partial
#' decomposition falls back to [base::eigen()] unless `force_engine = TRUE`.
#' The `"irlba"` option is accepted for compatibility but redirected because
#' an SVD is not a signed eigendecomposition.
#' @param A A finite symmetric square dense or sparse matrix.
#' @param d Positive number of eigencomponents to return.
#' @param only_values Return eigenvalues without eigenvectors.
#' @param scale_by Scale by `"none"`, matrix `"dimension"`, or entry
#'   `"sparsity"`.
#' @param use_laplacian Transform `A` with [graph_laplacian()] first.
#' @param engine Decomposition backend.
#' @param force_engine Stop instead of falling back when the selected partial
#'   backend is unavailable or fails.
#' @param order_by Order by signed `"value"` or absolute `"magnitude"`.
#' @param safe_d_multiplier Multiplier for extra partial components computed
#'   before selecting the requested components.
#' @param validate_inputs Validate matrix and option inputs.
#' @return A list with `values` and, unless `only_values`, `vectors`.
#' @examples eig_decomp(diag(c(3, 2, 1)), d = 2, engine = "base")
#' @seealso [singular_decomp()], [ase()]
#' @name eig_decomp
#' @export
eig_decomp <- function(
    A,
    d,
    only_values = FALSE,
    scale_by = c("none", "dimension", "sparsity"),
    use_laplacian = FALSE,
    engine = c("rspectra", "base", "irlba"),
    force_engine = FALSE,
    order_by = c("magnitude", "value"),
    safe_d_multiplier = 1,
    validate_inputs = TRUE) {
  if (length(validate_inputs) != 1L || !is.logical(validate_inputs) ||
      is.na(validate_inputs)) {
    stop("validate_inputs must be TRUE or FALSE.", call. = FALSE)
  }
  if (validate_inputs) {
    if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
      stop("A must be a square matrix-like object.", call. = FALSE)
    }
    if (nrow(A) < 1L) {
      stop("A must have positive dimensions.", call. = FALSE)
    }
    if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
      stop("A must contain only finite numeric values.", call. = FALSE)
    }
    if (!is_symmetric_matrix(A)) {
      stop("A must be symmetric.", call. = FALSE)
    }
  }
  if (length(d) != 1L ||
      !is.numeric(d) ||
      is.na(d) ||
      !is.finite(d) ||
      d < 1 ||
      d != floor(d)) {
    stop("d must be one positive integer.", call. = FALSE)
  }
  if (length(safe_d_multiplier) != 1L ||
      !is.numeric(safe_d_multiplier) ||
      is.na(safe_d_multiplier) ||
      !is.finite(safe_d_multiplier) ||
      safe_d_multiplier < 1) {
    stop("safe_d_multiplier must be one finite number at least 1.",
         call. = FALSE)
  }

  logical_inputs <- list(
    only_values = only_values,
    use_laplacian = use_laplacian,
    force_engine = force_engine
  )
  invalid_logical <- vapply(
    logical_inputs,
    function(value) {
      length(value) != 1L || !is.logical(value) || is.na(value)
    },
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_inputs)[invalid_logical], collapse = ", "),
      " must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  scale_by <- match.arg(scale_by)
  engine <- match.arg(engine)
  order_by <- match.arg(order_by)
  n <- nrow(A)
  d_requested <- min(as.integer(d), n)
  d_compute <- min(
    ceiling(safe_d_multiplier * d_requested),
    n
  )

  if (scale_by == "sparsity") {
    scale_factor <- n * matrix_density(A)
    if (!is.finite(scale_factor) ||
        abs(scale_factor) < .Machine$double.eps) {
      stop("Cannot scale by sparsity because matrix density is zero.",
           call. = FALSE)
    }
    A <- A / scale_factor
  } else if (scale_by == "dimension") {
    A <- A / n
  }
  if (use_laplacian) {
    A <- graph_laplacian(A)
  }

  base_eig <- function() {
    eig <- eigen(as.matrix(A), symmetric = TRUE, only.values = only_values)
    order_id <- if (order_by == "magnitude") {
      order(abs(eig$values), decreasing = TRUE)[seq_len(d_requested)]
    } else {
      order(eig$values, decreasing = TRUE)[seq_len(d_requested)]
    }
    out <- list(values = eig$values[order_id])
    if (!only_values) {
      out$vectors <- eig$vectors[, order_id, drop = FALSE]
    }
    out
  }

  if (engine == "irlba") {
    if (force_engine) {
      stop("irlba computes an SVD and is not a valid signed eigensolver.",
           call. = FALSE)
    }
    engine <- if (requireNamespace("RSpectra", quietly = TRUE)) {
      message("irlba is not a signed eigensolver; using rspectra.")
      "rspectra"
    } else {
      message("irlba is not a signed eigensolver; using base.")
      "base"
    }
  }
  if (n <= 200L && engine != "base" && !force_engine) {
    message("n <= 200; using engine = 'base'.")
    engine <- "base"
  }
  if (engine == "rspectra" && !requireNamespace("RSpectra", quietly = TRUE)) {
    if (force_engine) {
      stop("The RSpectra package is unavailable.", call. = FALSE)
    }
    message("RSpectra is unavailable; using engine = 'base'.")
    engine <- "base"
  }
  if (engine != "base" && d_compute >= n) {
    if (d_requested >= n) {
      message("A full eigendecomposition is required; using engine = 'base'.")
      engine <- "base"
    } else {
      d_compute <- n - 1L
    }
  }
  if (engine == "base") {
    return(base_eig())
  }

  tryCatch({
    eig <- RSpectra::eigs_sym(
      A,
      k = d_compute,
      which = if (order_by == "magnitude") "LM" else "LA",
      opts = list(retvec = !only_values)
    )
    order_id <- if (order_by == "magnitude") {
      order(abs(eig$values), decreasing = TRUE)[seq_len(d_requested)]
    } else {
      order(eig$values, decreasing = TRUE)[seq_len(d_requested)]
    }
    out <- list(values = eig$values[order_id])
    if (!only_values) {
      out$vectors <- eig$vectors[, order_id, drop = FALSE]
    }
    out
  }, error = function(e) {
    if (force_engine) {
      stop("RSpectra eigendecomposition failed: ", conditionMessage(e),
           call. = FALSE)
    }
    message("RSpectra eigendecomposition failed; using base: ",
            conditionMessage(e))
    base_eig()
  })
}

# Compute selected singular values and optional left/right singular vectors.
#
# Partial engines fall back to base::svd when unavailable or unsuccessful
# unless force_engine is TRUE. Singular values are non-negative, so the value
# and magnitude orderings are equivalent and both are retained for API parity.
#' Compute a singular-value decomposition
#'
#' Computes selected singular values and optional left and right vectors.
#' Partial backends fall back to [base::svd()] unless `force_engine = TRUE`.
#' @param A A finite dense or sparse matrix.
#' @param d Positive number of singular components to return.
#' @param only_values Return singular values without singular vectors.
#' @param nu,nv Numbers of left and right singular vectors to return.
#' @param scale_by Scale by `"none"`, matrix `"dimension"`, or entry
#'   `"sparsity"`.
#' @param use_laplacian Transform `A` with [graph_laplacian()] first.
#' @param engine Decomposition backend: `"rspectra"`, `"base"`, or `"irlba"`.
#' @param force_engine Stop instead of falling back when the selected partial
#'   backend is unavailable or fails.
#' @param order_by Component ordering convention; singular values make
#'   `"value"` and `"magnitude"` equivalent.
#' @param safe_d_multiplier Multiplier for extra partial components computed
#'   before selecting the requested components.
#' @param validate_inputs Validate matrix and option inputs.
#' @return A list with `values` and requested `u` and `v` matrices.
#' @examples singular_decomp(diag(c(3, 2, 1)), d = 2, engine = "base")
#' @seealso [eig_decomp()], [truncated_svd_reconstruct()]
#' @name singular_decomp
#' @export
singular_decomp <- function(
    A,
    d,
    only_values = FALSE,
    nu = if (only_values) 0L else d,
    nv = if (only_values) 0L else d,
    scale_by = c("none", "dimension", "sparsity"),
    use_laplacian = FALSE,
    engine = c("rspectra", "base", "irlba"),
    force_engine = FALSE,
    order_by = c("value", "magnitude"),
    safe_d_multiplier = 1,
    validate_inputs = TRUE) {
  if (length(validate_inputs) != 1L || !is.logical(validate_inputs) ||
      is.na(validate_inputs)) {
    stop("validate_inputs must be TRUE or FALSE.", call. = FALSE)
  }
  if (validate_inputs) {
    if (is.null(dim(A)) || length(dim(A)) != 2L || any(dim(A) < 1L)) {
      stop("A must be a matrix-like object with positive dimensions.",
           call. = FALSE)
    }
    if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
      stop("A must contain only finite numeric values.", call. = FALSE)
    }
  }
  logical_inputs <- list(
    only_values = only_values,
    use_laplacian = use_laplacian,
    force_engine = force_engine
  )
  invalid_logical <- vapply(
    logical_inputs,
    function(value) {
      length(value) != 1L || !is.logical(value) || is.na(value)
    },
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_inputs)[invalid_logical], collapse = ", "),
      " must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  min_dim <- min(dim(A))
  integer_inputs <- list(d = d, nu = nu, nv = nv)
  invalid_integer <- vapply(
    integer_inputs,
    function(value) {
      length(value) != 1L ||
        !is.numeric(value) ||
        is.na(value) ||
        !is.finite(value) ||
        value < 0 ||
        value != floor(value)
    },
    logical(1)
  )
  invalid_integer["d"] <- invalid_integer["d"] || d < 1
  if (any(invalid_integer)) {
    stop(
      paste(names(integer_inputs)[invalid_integer], collapse = ", "),
      " must be valid integer dimensions (with d positive).",
      call. = FALSE
    )
  }
  if (length(safe_d_multiplier) != 1L ||
      !is.numeric(safe_d_multiplier) ||
      is.na(safe_d_multiplier) ||
      !is.finite(safe_d_multiplier) ||
      safe_d_multiplier < 1) {
    stop("safe_d_multiplier must be one finite number at least 1.",
         call. = FALSE)
  }
  scale_by <- match.arg(scale_by)
  engine <- match.arg(engine)
  order_by <- match.arg(order_by)
  n <- nrow(A)
  p <- ncol(A)
  d_requested <- min(as.integer(d), min_dim)
  d_compute <- min(ceiling(safe_d_multiplier * d_requested), min_dim)
  nu_requested <- if (only_values) 0L else min(as.integer(nu), d_requested, n)
  nv_requested <- if (only_values) 0L else min(as.integer(nv), d_requested, p)

  if (scale_by == "sparsity") {
    scale_factor <- sqrt(as.double(n) * p) * matrix_density(A)
    if (!is.finite(scale_factor) ||
        abs(scale_factor) < .Machine$double.eps) {
      stop("Cannot scale by sparsity because matrix density is zero.",
           call. = FALSE)
    }
    A <- A / scale_factor
  } else if (scale_by == "dimension") {
    A <- A / sqrt(as.double(n) * p)
  }
  if (use_laplacian) {
    A <- graph_laplacian(A)
  }

  base_svd <- function() {
    sv <- svd(as.matrix(A), nu = nu_requested, nv = nv_requested)
    order_id <- order(sv$d, decreasing = TRUE)[seq_len(d_requested)]
    out <- list(values = sv$d[order_id])
    if (!only_values && nu_requested > 0L) {
      out$u <- sv$u[, order_id[seq_len(nu_requested)], drop = FALSE]
    }
    if (!only_values && nv_requested > 0L) {
      out$v <- sv$v[, order_id[seq_len(nv_requested)], drop = FALSE]
    }
    out
  }

  package_name <- if (engine == "rspectra") "RSpectra" else engine
  if (min_dim <= 200L && engine != "base" && !force_engine) {
    message("min(dim(A)) <= 200; using engine = 'base'.")
    engine <- "base"
  } else if (engine != "base" &&
             !requireNamespace(package_name, quietly = TRUE)) {
    if (force_engine) {
      stop("The ", package_name, " package is unavailable.", call. = FALSE)
    }
    message(package_name, " is unavailable; using engine = 'base'.")
    engine <- "base"
  }
  if (engine != "base" && d_compute >= min_dim) {
    if (d_requested >= min_dim) {
      message("A full singular decomposition is required; using engine = 'base'.")
      engine <- "base"
    } else {
      d_compute <- min_dim - 1L
    }
  }
  if (engine == "base") {
    return(base_svd())
  }

  partial_svd <- function() {
    need_u <- !only_values && nu_requested > 0L
    need_v <- !only_values && nv_requested > 0L
    compute_nu <- if (need_u) d_compute else 0L
    compute_nv <- if (need_v || !need_u) d_compute else 0L
    if (engine == "irlba") {
      irlba::irlba(A, nu = compute_nu, nv = compute_nv)
    } else {
      RSpectra::svds(A, k = d_compute, nu = compute_nu, nv = compute_nv)
    }
  }

  tryCatch({
    sv <- partial_svd()
    order_id <- order(sv$d, decreasing = TRUE)[seq_len(d_requested)]
    out <- list(values = sv$d[order_id])
    if (!only_values && nu_requested > 0L) {
      out$u <- sv$u[, order_id[seq_len(nu_requested)], drop = FALSE]
    }
    if (!only_values && nv_requested > 0L) {
      out$v <- sv$v[, order_id[seq_len(nv_requested)], drop = FALSE]
    }
    out
  }, error = function(e) {
    if (force_engine) {
      stop(engine, " singular decomposition failed: ", conditionMessage(e),
           call. = FALSE)
    }
    message(engine, " singular decomposition failed; using base: ",
            conditionMessage(e))
    base_svd()
  })
}

# Reconstruct and clip the rank-d truncated singular value decomposition of A.
#
# Singular values scale the columns of the left singular-vector matrix using a
# vectorized sweep. tcrossprod then forms U_d diag(s_d) t(V_d) without creating
# a diagonal matrix or using the general matrix-multiplication operator. The
# default [0, 1] bounds produce a probability-matrix estimate. Use infinite
# bounds to obtain the unmodified truncated-SVD reconstruction.
#' Reconstruct a truncated singular-value approximation
#'
#' Reconstructs `U_d diag(s_d) V_d'` without materializing the diagonal matrix,
#' then clips entries to the requested interval. The default bounds produce a
#' probability-matrix estimate; infinite bounds leave the reconstruction
#' unmodified.
#' @param A A finite dense or sparse matrix.
#' @param d Positive reconstruction rank.
#' @param lower_clip,upper_clip Inclusive clipping bounds.
#' @param engine Singular-decomposition backend.
#' @param force_engine Stop instead of falling back when the backend fails.
#' @param safe_d_multiplier Multiplier for extra partial components.
#' @return A numeric rank-`d` matrix approximation.
#' @examples truncated_svd_reconstruct(diag(3), d = 1, engine = "base")
#' @seealso [singular_decomp()], [usvt()]
#' @name truncated_svd_reconstruct
#' @export
truncated_svd_reconstruct <- function(
    A,
    d,
    lower_clip = 0,
    upper_clip = 1,
    engine = c("rspectra", "base", "irlba"),
    force_engine = FALSE,
    safe_d_multiplier = 1) {
  engine <- match.arg(engine)
  decomposition <- singular_decomp(
    A = A,
    d = d,
    only_values = FALSE,
    nu = d,
    nv = d,
    scale_by = "none",
    use_laplacian = FALSE,
    engine = engine,
    force_engine = force_engine,
    order_by = "value",
    safe_d_multiplier = safe_d_multiplier
  )

  scaled_u <- sweep(
    decomposition$u,
    2L,
    decomposition$values,
    `*`
  )
  estimate <- tcrossprod(scaled_u, decomposition$v)
  dimnames(estimate) <- dimnames(A)
  clip_values(
    estimate,
    lower_clip = lower_clip,
    upper_clip = upper_clip
  )
}

# Estimate a probability matrix using the pasted USVT rank rule.
#
# This historical variant sets d = ceiling(n^(1/3)), computes the corresponding
# truncated SVD, and clips every reconstructed entry to [0, 1]. Its public name
# is shortened to `usvt` under the project naming rules.
#' Estimate a matrix with universal singular-value thresholding
#'
#' Applies the historical pasted USVT rule `d = ceiling(n^(1/3))`, reconstructs
#' that truncated SVD, and clips every entry to the supplied interval.
#' @param A A finite nonempty square dense or sparse matrix.
#' @param lower_clip,upper_clip Inclusive clipping bounds.
#' @param engine Singular-decomposition backend.
#' @param force_engine Stop instead of falling back when the backend fails.
#' @param safe_d_multiplier Multiplier for extra partial components.
#' @return A dense thresholded matrix estimate.
#' @examples usvt(diag(3), engine = "base")
#' @seealso [singular_decomp()], [ase()]
#' @name usvt
#' @export
usvt <- function(
    A,
    lower_clip = 0,
    upper_clip = 1,
    engine = c("irlba", "rspectra", "base"),
    force_engine = FALSE,
    safe_d_multiplier = 1) {
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
    stop("A must be a square matrix-like object.", call. = FALSE)
  }
  if (nrow(A) < 1L) {
    stop("A must have positive dimensions.", call. = FALSE)
  }

  engine <- match.arg(engine)
  n <- nrow(A)
  d_usvt <- min(ceiling(n^(1 / 3)), n)
  P_hat <- truncated_svd_reconstruct(
    A = A,
    d = d_usvt,
    lower_clip = lower_clip,
    upper_clip = upper_clip,
    engine = engine,
    force_engine = force_engine,
    safe_d_multiplier = safe_d_multiplier
  )
  P_hat
}
