# ============================================================
# Structure-preserving calibration of binary-outcome DGMs
# ============================================================
#
# DISCLAIMER
# -------
# AI assistance: Microsoft Copilot was used to assist with annotating 
# the R code and drafting the accompanying documentation. 
# The underlying methodology and implementation were developed by the 
# authors, who reviewed and take responsibility for the final materials.
#
#
#
# PURPOSE
# -------
# This script defines the functions used to construct and inspect a
# structure-preserving calibration path for binary-outcome data-generating
# mechanisms (DGMs). Along this path, the directions of the arm-specific
# linear predictors are held fixed, while a common scaling parameter
# (lambda) controls their overall strength. Arm-specific intercepts are
# recalibrated so that the marginal outcome risks remain equal to the
# user-specified values mu0 and mu1.
#
# The resulting path can be indexed by the relative marginal prognostic
# strength (MPS) in the control arm. At each point, the script evaluates:
#   * control-arm and treatment-arm relative MPS;
#   * the corresponding scaling parameter and calibrated intercepts; and
#   * precision-gain potential (PGP) for the risk difference (RD), risk
#     ratio (RR), marginal odds ratio (OR), and optional additional smooth
#     marginal estimands.
#
# USAGE
# -----
# Save this file, then make all functions available in an R session with:
#
#   source("precision_gain_functions.R")
#
#
# NOTATION AND INPUT CONVENTIONS
# ------------------------------
# h       : N x p matrix of evaluated covariate features for a reference
#           sample. Rows represent observations; columns represent features.
# v_hom   : length-p coefficient-direction vector shared by both arms.
# v_het   : length-p treatment-arm deviation from v_hom.
# eta0    : control-arm score h %*% v_hom.
# eta1    : treatment-arm score h %*% (v_hom + v_het).
# mu0     : target marginal outcome risk under control.
# mu1     : target marginal outcome risk under treatment.
# pi      : treatment-allocation probability P(A = 1).
# lambda  : nonnegative common scale applied to eta0 and eta1.
# alpha0,
# alpha1  : arm-specific intercepts calibrated to preserve mu0 and mu1.
#
# All empirical expectations and variances below are computed over the
# supplied reference sample.


# ------------------------------------------------------------
# 1. Relative MPS at a fixed value of lambda
# ------------------------------------------------------------
# For a supplied score eta and marginal risk mu, calibrate the intercept
# alpha so that the empirical mean of
#
#   p(X) = expit(alpha + lambda * eta)
#
# equals mu. The relative MPS is then
#
#   mean[{p(X) - mu}^2] / {mu * (1 - mu)}.
#
# Args:
#   lambda: Nonnegative scale parameter for the score eta.
#   eta:    Numeric vector containing the score in one treatment arm.
#   mu:     Target marginal outcome risk; must lie strictly between 0 and 1.
#
# Returns:
#   A list with the calibrated intercept alpha, the vector of conditional
#   risks p, and the resulting relative MPS.
mps <- function(lambda, eta, mu) {
  # The wide bracket is intended to cover all practically relevant
  # intercepts. uniroot() finds the intercept enforcing mean(p) = mu.
  alpha <- uniroot(
    function(alpha) mean(plogis(alpha + lambda * eta)) - mu,
    c(-1000, 1000)
  )$root

  p <- plogis(alpha + lambda * eta)

  list(
    alpha = alpha,
    p = p,
    mps = mean((p - mu)^2) / (mu * (1 - mu))
  )
}


# ------------------------------------------------------------
# 2. Empirical limiting relative MPS
# ------------------------------------------------------------
# Approximate the largest relative MPS reached as lambda tends to infinity,
# while the intercept continues to preserve the marginal risk mu.
#
# In that limit, observations with the largest score values receive risk 1,
# observations with smaller scores receive risk 0, and observations tied at
# the threshold may receive a common fractional risk q so that the empirical
# mean remains exactly mu. Explicit tie handling is important because eta may
# contain repeated values, particularly with discrete covariates.
#
# Args:
#   eta: Numeric score vector for one treatment arm.
#   mu:  Target marginal outcome risk.
#
# Returns:
#   A single numeric value: the empirical limiting relative MPS.
limiting_mps <- function(eta, mu) {
  # Sort scores from largest to smallest and collect ties.
  groups <- rle(sort(eta, decreasing = TRUE))
  cumulative <- cumsum(groups$lengths)

  # Number of outcome-risk units required to obtain empirical mean mu.
  target <- mu * length(eta)

  # First tied score group at which the required total is reached.
  threshold <- which(cumulative >= target)[1]
  full <- if (threshold == 1) 0 else cumulative[threshold - 1]

  # Fractional risk assigned to observations in the threshold group.
  q <- (target - full) / groups$lengths[threshold]

  # E[p(X)^2] - mu^2 for the limiting step function.
  variance_limit <-
    (full + groups$lengths[threshold] * q^2) / length(eta) - mu^2

  variance_limit / (mu * (1 - mu))
}


# ------------------------------------------------------------
# 3. Calibrate lambda to a target relative MPS
# ------------------------------------------------------------
# Find the nonnegative value of lambda that produces a requested relative
# MPS. Every evaluation of the outer root-finding objective calls MPS(),
# which internally recalibrates the intercept to preserve mu.
#
# Args:
#   target_mps: Requested relative MPS.
#   eta:        Numeric score vector for the arm used for calibration.
#   mu:         Target marginal risk in that arm.
#
# Returns:
#   A list containing lambda and the corresponding alpha, p, and mps.
calibrate_lambda <- function(target_mps, eta, mu) {
  if (target_mps == 0) {
    # At lambda = 0, the conditional risk is constant and MPS is zero.
    lambda <- 0
  } else {
    objective <- function(lambda) {
      mps(lambda, eta, mu)$mps - target_mps
    }

    # Expand the upper endpoint until it brackets the requested MPS.
    upper <- 1
    while (objective(upper) < 0) {
      upper <- 2 * upper
    }

    lambda <- uniroot(objective, c(0, upper))$root
  }

  c(list(lambda = lambda), mps(lambda, eta, mu))
}


# ------------------------------------------------------------
# 4. Define marginal estimands through their gradients
# ------------------------------------------------------------
# PGP depends on the gradient of a smooth marginal estimand f(mu1, mu0).
# This helper always supplies gradients for RD, RR, and the marginal OR.
# Users may add other differentiable estimands by supplying named list
# entries with a gradient function and, optionally, a display label.
#
# Gradient convention:
#   c(partial derivative with respect to mu1,
#     partial derivative with respect to mu0)
#
# Args:
#   mu0, mu1:            Control- and treatment-arm marginal risks.
#   additional_estimands: NULL, or a named list. Each entry must contain
#                         gradient = function(mu1, mu0) ..., returning two
#                         finite values. An optional label may be supplied.
#
# Returns:
#   A named list. Each entry contains a display label and a length-two
#   gradient vector evaluated at the supplied marginal risks.
estimands <- function(mu0, mu1, additional_estimands = NULL) {
  marginal_or <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))

  out <- list(
    RD = list(
      label = "Risk difference",
      gradient = c(1, -1)
    ),
    RR = list(
      label = "Risk ratio",
      gradient = c(1 / mu0, -mu1 / mu0^2)
    ),
    OR = list(
      label = "Odds ratio",
      gradient = marginal_or * c(
        1 / (mu1 * (1 - mu1)),
        -1 / (mu0 * (1 - mu0))
      )
    )
  )

  if (!is.null(additional_estimands)) {
    for (name in names(additional_estimands)) {
      specification <- additional_estimands[[name]]
      gradient <- specification$gradient(mu1, mu0)

      if (length(gradient) != 2 || any(!is.finite(gradient))) {
        stop(
          "Each additional estimand must return a finite gradient of length two."
        )
      }

      out[[name]] <- list(
        label = if (is.null(specification$label)) name else specification$label,
        gradient = gradient
      )
    }
  }

  out
}


# ------------------------------------------------------------
# 5. Evaluate one DGM on the calibration path
# ------------------------------------------------------------
# The path is indexed by control-arm relative MPS. First, lambda is calibrated
# using eta0 and mu0. The same lambda is then applied to eta1, with a separate
# treatment-arm intercept calibrated to preserve mu1. This common scaling is
# what preserves the specified coefficient directions across the path.
#
# Args:
#   target_mps0:          Target control-arm relative MPS.
#   eta0, eta1:           Control- and treatment-arm score vectors.
#   mu0, mu1:             Target marginal risks.
#   pi:                   Treatment-allocation probability.
#   additional_estimands: Optional estimand specifications passed to
#                         estimands().
#
# Returns:
#   A list containing both relative MPS values, lambda, both intercepts,
#   estimand-specific PGP values, and the evaluated estimand definitions.
evaluate_dgm <- function(
    target_mps0,
    eta0,
    eta1,
    mu0,
    mu1,
    pi = 0.5,
    additional_estimands = NULL) {

  # Calibrate lambda using the control arm, then apply the same lambda in
  # the treatment arm and recalibrate only its intercept.
  fit0 <- calibrate_lambda(target_mps0, eta0, mu0)
  fit1 <- mps(fit0$lambda, eta1, mu1)

  # Empirical arm-specific prognostic variances and their cross-moment.
  S0 <- mean((fit0$p - mu0)^2)
  S1 <- mean((fit1$p - mu1)^2)
  C01 <- mean((fit0$p - mu0) * (fit1$p - mu1))

  estimands_list <- estimands(mu0, mu1, additional_estimands)

  # Evaluate PGP separately for each estimand. The numerator is the
  # first-order variance reduction, and the denominator is the asymptotic
  # variance under unadjusted estimation.
  pgp <- vapply(estimands_list, function(estimand) {
    f1 <- estimand$gradient[1]
    f0 <- estimand$gradient[2]

    gain <-
      (1 - pi) / pi * f1^2 * S1 +
      pi / (1 - pi) * f0^2 * S0 -
      2 * f1 * f0 * C01

    variance_unadjusted <-
      f1^2 * mu1 * (1 - mu1) / pi +
      f0^2 * mu0 * (1 - mu0) / (1 - pi)

    gain / variance_unadjusted
  }, numeric(1))

  list(
    control_mps = fit0$mps,
    treatment_mps = fit1$mps,
    lambda = fit0$lambda,
    alpha0 = fit0$alpha,
    alpha1 = fit1$alpha,
    pgp = pgp,
    estimands = estimands_list
  )
}


# ------------------------------------------------------------
# 6. Construct a finite display grid over the calibration path
# ------------------------------------------------------------
# Build arm-specific scores from the reference-sample feature matrix and the
# coefficient vectors. The finite grid runs from zero to a fraction
# t of the empirical limiting control-arm MPS. It is intended for display and
# exploration; select_dgm() below performs exact recalibration at a chosen MPS
# and therefore does not depend on the grid resolution K.
#
# Args:
#   h:                     N x p feature matrix.
#   v_hom:                 Length-p shared coefficient-direction vector.
#   v_het:                 Length-p treatment-arm deviation vector.
#   mu0, mu1:              Target marginal risks.
#   pi:                    Treatment-allocation probability.
#   t:                     Fraction of the limiting control-arm MPS included
#                          in the displayed grid; typically below 1.
#   K:                     Number of grid points.
#   additional_estimands:  Optional estimand specifications.
#
# Returns:
#   A "basis" list. Its path component is a data frame with one row per grid
#   point and columns for both MPS values, lambda, both intercepts, and all
#   requested PGP values. The remaining components retain the information
#   needed for plotting and exact recalibration.
calibration_path <- function(
    h,
    v_hom,
    v_het,
    mu0,
    mu1,
    pi = 0.5,
    t = 0.95,
    K = 100,
    additional_estimands = NULL) {

  # Arm-specific unscaled scores. drop() converts one-column matrix products
  # into ordinary numeric vectors.
  eta0 <- drop(h %*% v_hom)
  eta1 <- drop(h %*% (v_hom + v_het))

  if (var(eta0) == 0) {
    stop("The control-arm score must be nonconstant.")
  }

  # The grid deliberately stops at t times the empirical limiting value to
  # avoid requiring an infinite lambda at the endpoint.
  s0_infinity <- limiting_mps(eta0, mu0)
  s0_grid <- seq(0, t * s0_infinity, length.out = K)

  # Evaluate every grid point while keeping all non-indexing inputs fixed.
  fits <- lapply(
    s0_grid,
    evaluate_dgm,
    eta0 = eta0,
    eta1 = eta1,
    mu0 = mu0,
    mu1 = mu1,
    pi = pi,
    additional_estimands = additional_estimands
  )

  # Collect scalar calibration outputs.
  path <- data.frame(
    control_mps = vapply(fits, `[[`, numeric(1), "control_mps"),
    treatment_mps = vapply(fits, `[[`, numeric(1), "treatment_mps"),
    lambda = vapply(fits, `[[`, numeric(1), "lambda"),
    alpha0 = vapply(fits, `[[`, numeric(1), "alpha0"),
    alpha1 = vapply(fits, `[[`, numeric(1), "alpha1")
  )

  # Add one PGP column per estimand.
  pgp <- t(vapply(
    fits,
    `[[`,
    numeric(length(fits[[1]]$pgp)),
    "pgp"
  ))
  colnames(pgp) <- names(fits[[1]]$pgp)
  path <- cbind(path, as.data.frame(pgp, check.names = FALSE))

  list(
    path = path,
    eta0 = eta0,
    eta1 = eta1,
    v_hom = v_hom,
    v_het = v_het,
    mu0 = mu0,
    mu1 = mu1,
    pi = pi,
    t = t,
    s0_infinity = s0_infinity,
    estimands = fits[[1]]$estimands,
    additional_estimands = additional_estimands
  )
}


# ------------------------------------------------------------
# 7. Recalibrate exactly at a user-selected control-arm MPS
# ------------------------------------------------------------
# Select a DGM at any admissible control-arm relative MPS. The function calls
# evaluate_dgm() directly, rather than interpolating between rows of basis$path.
# Consequently, its result is independent of the display-grid resolution K.
#
# Args:
#   basis:       Object returned by calibration_path().
#   target_mps0: Desired control-arm relative MPS within the displayed range.
#
# Returns:
#   The output of evaluate_dgm(), augmented with the final shared and
#   heterogeneous coefficient vectors beta_hom and beta_het.
select_dgm <- function(basis, target_mps0) {
  if (target_mps0 < 0 || target_mps0 > basis$t * basis$s0_infinity) {
    stop("target_mps0 lies outside the grid range.")
  }

  fit <- evaluate_dgm(
    target_mps0,
    basis$eta0,
    basis$eta1,
    basis$mu0,
    basis$mu1,
    basis$pi,
    basis$additional_estimands
  )

  # Lambda converts the coefficient directions into the coefficients of the
  # selected DGM. The treatment-arm coefficient vector is their sum.
  fit$beta_hom <- fit$lambda * basis$v_hom
  fit$beta_het <- fit$lambda * basis$v_het

  fit
}


# ------------------------------------------------------------
# 8. Plot the calibration path
# ------------------------------------------------------------
# Produce a one-panel base-R display of treatment-arm relative MPS and the
# requested PGP curves against control-arm relative MPS. All displayed values
# are percentages. The diagonal reference line marks equal MPS in both arms.
#
# Args:
#   basis:          Object returned by calibration_path().
#   file:           Optional PDF filename. If NULL, draw on the active device.
#   show_estimands: Character vector naming PGP curves to display.
#   selected_mps:   Optional control-arm relative MPS to mark in the plot.
#   tick_by:        Spacing of labeled ticks, in percentage points.
#   grid_by:        Spacing of background grid lines, in percentage points.
#   estimand_lwd:   Line width for PGP curves.
#   mps_lwd:        Line width for the treatment-arm MPS curve.
#
# Returns:
#   NULL invisibly. The function's purpose is its plotting side effect.
plot_calibration_path <- function(
    basis,
    file = NULL,
    show_estimands = names(basis$estimands),
    selected_mps = NULL,
    tick_by = 10,
    grid_by = 5,
    estimand_lwd = 1.25,
    mps_lwd = 2.25) {

  # If a filename is supplied, direct the plot to a PDF and guarantee that
  # the device is closed when the function exits, including after an error.
  if (!is.null(file)) {
    pdf(file, width = 8, height = 6, family = "sans")
    on.exit(dev.off(), add = TRUE)
  }

  path <- basis$path

  # Fail early if the user requests an estimand not contained in basis.
  unknown_estimands <- setdiff(show_estimands, names(basis$estimands))

  if (length(unknown_estimands) > 0) {
    stop(
      "Unknown estimand(s): ",
      paste(unknown_estimands, collapse = ", "),
      ". Available estimands are: ",
      paste(names(basis$estimands), collapse = ", ")
    )
  }

  if (length(show_estimands) == 0) {
    stop("At least one estimand must be selected.")
  }

  # Convert all displayed quantities from proportions to percentages.
  x <- 100 * path$control_mps
  treatment_mps <- 100 * path$treatment_mps
  pgp <- 100 * as.matrix(path[, show_estimands, drop = FALSE])

  labels <- vapply(
    basis$estimands[show_estimands],
    `[[`,
    character(1),
    "label"
  )

  estimand_colors <- hcl.colors(length(show_estimands), "Dark 3")

  # Shared limits make the identity line and all curves directly comparable.
  x_limit <- ceiling(max(x) / tick_by) * tick_by
  y_limit <- ceiling(max(c(x, treatment_mps, pgp)) / tick_by) * tick_by
  x_limit <- max(x_limit, tick_by)
  y_limit <- max(y_limit, tick_by)

  x_ticks <- seq(from = 0, to = x_limit, by = tick_by)
  y_ticks <- seq(from = 0, to = y_limit, by = tick_by)
  x_grid <- seq(from = 0, to = x_limit, by = grid_by)
  y_grid <- seq(from = 0, to = y_limit, by = grid_by)

  # Preserve and restore the user's graphics settings.
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(4.5, 4.7, 1.0, 1.0), las = 1)

  # Initialize an empty plotting region so the background grid can be drawn
  # before the data curves.
  plot(
    x,
    treatment_mps,
    type = "n",
    xlab = "Control-arm relative MPS (%)",
    ylab = "Relative marginal property (%)",
    xlim = c(0, x_limit),
    ylim = c(0, y_limit),
    xaxt = "n",
    yaxt = "n"
  )

  axis(side = 1, at = x_ticks, labels = x_ticks)
  axis(side = 2, at = y_ticks, labels = y_ticks)

  abline(h = y_grid, col = "grey90", lty = 1, lwd = 0.6)
  abline(v = x_grid, col = "grey90", lty = 1, lwd = 0.6)

  # Identity line: control-arm MPS equals treatment-arm MPS.
  abline(a = 0, b = 1, lty = 2, lwd = 1, col = "grey45")

  # Treatment-arm MPS and estimand-specific PGP curves.
  lines(x, treatment_mps, lwd = mps_lwd, col = "#31688E")
  matlines(x, pgp, lty = 1, lwd = estimand_lwd, col = estimand_colors)

  # Optionally mark a selected control-arm MPS. Values on the displayed
  # curves are obtained by linear interpolation for plotting only. Exact DGM
  # recalibration is performed by select_dgm().
  if (!is.null(selected_mps)) {
    selected_x <- 100 * selected_mps

    if (selected_x < 0 || selected_x > max(x)) {
      warning("selected_mps lies outside the displayed calibration range.")
    } else {
      abline(v = selected_x, lty = 3, lwd = 1.25, col = "black")

      points(
        selected_x,
        approx(x = x, y = treatment_mps, xout = selected_x)$y,
        pch = 19,
        cex = 0.8,
        col = "#31688E"
      )

      for (j in seq_along(show_estimands)) {
        points(
          selected_x,
          approx(x = x, y = pgp[, j], xout = selected_x)$y,
          pch = 19,
          cex = 0.7,
          col = estimand_colors[j]
        )
      }
    }
  }

  legend_labels <- c("Treatment-arm MPS", labels, "Identity line")
  legend_colors <- c("#31688E", estimand_colors, "grey45")
  legend_lty <- c(1, rep(1, length(show_estimands)), 2)
  legend_lwd <- c(mps_lwd, rep(estimand_lwd, length(show_estimands)), 1)

  legend(
    "topleft",
    legend = legend_labels,
    col = legend_colors,
    lty = legend_lty,
    lwd = legend_lwd,
    bty = "n"
  )

  invisible(NULL)
}
