#' Gaussian Mixture EM -- direct R port of GaussianMixture.m (Bowman CLUSTER)
#'
#' Fits a K-component univariate Gaussian mixture using the same EM algorithm
#' bundled with CharAnalysis v1.1 and v2.0 (Bowman CLUSTER implementation).
#' Replicates three key behaviours that distinguish it from \code{mclust}:
#'
#' \enumerate{
#'   \item \strong{First/last-point initialisation}: for K = 2, component
#'     means are seeded at the first and last elements of the input data
#'     vector (not at sorted min/max; the data are in time order).  This is
#'     the \code{initMixture()} behaviour from MATLAB.
#'   \item \strong{Loose convergence criterion}: EM stops when the per-step
#'     log-likelihood gain falls to or below
#'     \eqn{\epsilon = 0.01 \times L_c \times \log N} where \eqn{L_c = 3}
#'     for univariate data.  This is much looser than the \code{mclust}
#'     default (\eqn{10^{-5}}) and causes MATLAB to freeze closer to the
#'     initial configuration.
#'   \item \strong{Variance regularisation}: a small floor
#'     \eqn{R_{\min} = \bar{\sigma}^2 / 10^5} is added to each component
#'     variance after every M-step, preventing degenerate (zero-variance)
#'     solutions.
#' }
#'
#' Because CharAnalysis always calls \code{GaussianMixture(X, 2, 2, false)}
#' (i.e.\ \code{finalK = initK = 2}), the \code{MDLReduceOrder()} path is
#' never exercised and is therefore not implemented here.
#'
#' @param x  Numeric vector of observations (the C_peak values in the window).
#' @param k  Integer number of components (default 2).
#'
#' @return Named list:
#'   \describe{
#'     \item{mu}{Numeric vector length \code{k}: component means sorted
#'       ascending.}
#'     \item{sigma}{Numeric vector length \code{k}: component standard
#'       deviations (= sqrt of fitted variance), sorted to match \code{mu}.}
#'     \item{prop}{Numeric vector length \code{k}: mixing proportions, sorted
#'       to match \code{mu}.}
#'     \item{loglik}{Scalar: log-likelihood at convergence.}
#'   }
#'
#' @seealso [char_thresh_local()], [char_thresh_global()]
gaussian_mixture_em <- function(x, k = 2L) {

  x <- as.numeric(x)
  n <- length(x)
  k <- as.integer(k)

  # ---- Convergence threshold -------------------------------------------------
  # MATLAB's EMIterate uses epsilon = 0.01 * Lc * log(N*M) where for
  # univariate data (M=1) Lc = 1 + M + 0.5*M*(M+1) = 3, giving
  # epsilon = 0.03 * log(N).  For a typical 50-sample window this is ~0.117,
  # which is deliberately loose (Bowman CLUSTER convention) and causes the EM
  # to stop after only 1-2 iterations.  A tight criterion (1e-6) runs many
  # more iterations and diverges from MATLAB's solution, producing different
  # threshold values and therefore different peak counts.  Using MATLAB's
  # criterion keeps the two implementations in close agreement.
  eps <- 0.01 * 3 * log(n)         # MATLAB loose criterion: 0.03 * log(N)

  # ---- initMixture -----------------------------------------------------------
  # Population variance: MATLAB computes R = (N-1)/N * cov(pixels)
  # In R, var(x) is the sample variance (/(N-1)), so:
  r_init <- (n - 1L) / n * stats::var(x)
  rmin   <- r_init / 1e5           # regularisation floor (Rmin in MATLAB)

  # Evenly-spaced index initialisation:
  #   k=1 -> index 1  (first data point)
  #   k=2 -> index N  (last  data point) because period = (N-1)/(K-1)
  #   k=j -> index floor((j-1)*period) + 1
  period   <- (n - 1L) / (k - 1L)
  init_idx <- vapply(seq_len(k),
                     function(j) floor((j - 1L) * period) + 1L,
                     integer(1L))
  mu_k  <- x[init_idx]             # initial means
  r_k   <- rep(r_init + rmin, k)  # initial variances (equal at start)
  pb_k  <- rep(1 / k, k)          # equal initial mixing proportions

  # ---- EMIterate -------------------------------------------------------------
  ll_new <- .gm_loglik(x, mu_k, r_k, pb_k)

  repeat {
    ll_old <- ll_new

    # E-step: posterior responsibilities [N x K]
    resp <- .gm_responsibilities(x, mu_k, r_k, pb_k)

    # M-step
    nk   <- colSums(resp)
    mu_k <- colSums(x * resp) / nk
    for (j in seq_len(k)) {
      r_k[j] <- sum((x - mu_k[j])^2 * resp[, j]) / nk[j] + rmin
    }
    pb_k <- nk / sum(nk)           # ClusterNormalize: renormalise proportions

    ll_new <- .gm_loglik(x, mu_k, r_k, pb_k)

    # MATLAB stopping rule: break when gain <= epsilon (including negative gain)
    if ((ll_new - ll_old) <= eps) break
  }

  # ---- Sort by mean (ascending) and return ----------------------------------
  ord <- order(mu_k)
  list(
    mu     = mu_k[ord],
    sigma  = sqrt(r_k[ord]),
    prop   = pb_k[ord],
    loglik = ll_new
  )
}

# ---------------------------------------------------------------------------
# Internal helpers  (not exported; prefixed with .gm_ to avoid name clashes)
# ---------------------------------------------------------------------------

# Total log-likelihood of a K-component 1D Gaussian mixture at data x.
# Uses the log-sum-exp trick for numerical stability (same as MATLAB EStep).
.gm_loglik <- function(x, mu, r, pb) {
  n <- length(x)
  k <- length(mu)

  # log-density of each component at each observation
  log_pk <- matrix(NA_real_, n, k)
  for (j in seq_len(k)) {
    log_pk[, j] <- -0.5 * (x - mu[j])^2 / r[j] -
                   0.5 * (log(2 * pi) + log(r[j]))
  }

  # Log-sum-exp: subtract per-row maximum before exponentiating
  llmax <- apply(log_pk, 1L, max)
  pk    <- exp(log_pk - llmax) * matrix(pb, n, k, byrow = TRUE)
  ss    <- rowSums(pk)
  sum(log(ss) + llmax)
}

# Posterior responsibilities p(z_n = k | x_n) -- [N x K] matrix.
.gm_responsibilities <- function(x, mu, r, pb) {
  n <- length(x)
  k <- length(mu)

  log_pk <- matrix(NA_real_, n, k)
  for (j in seq_len(k)) {
    log_pk[, j] <- -0.5 * (x - mu[j])^2 / r[j] -
                   0.5 * (log(2 * pi) + log(r[j]))
  }

  llmax <- apply(log_pk, 1L, max)
  pk    <- exp(log_pk - llmax) * matrix(pb, n, k, byrow = TRUE)
  ss    <- rowSums(pk)
  pk / ss
}
