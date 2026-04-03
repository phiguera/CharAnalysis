#' Locally-weighted linear regression matching MATLAB charLowess.m / smooth()
#'
#' Pure-R implementation of the shifted-window lowess algorithm from MATLAB's
#' \code{charLowess.m} v2.0, which itself replicates MATLAB's
#' \code{smooth(y, k, 'lowess')} from the Curve Fitting Toolbox.  Replaces
#' the previous wrapper around \code{\link[stats]{lowess}}, which used a
#' shrinking window at record boundaries and different tricubic weight
#' normalisation, causing systematic differences of ~1e-3 in \code{charBkg}.
#'
#' @param y     Numeric vector of values to smooth (NaN-free; bridge NaNs
#'   before calling and restore afterward, matching the MATLAB workflow).
#' @param x     Numeric vector of x-positions (same length as \code{y}).
#'   If \code{NULL}, \code{seq_along(y)} is used (index coordinates).
#' @param span  Smoothing span.
#'   \itemize{
#'     \item \code{span >= 1}: number of data points.
#'     \item \code{span < 1}: fraction of data length.
#'   }
#'   Converted to integer \code{k = max(3, round(span * n))} or
#'   \code{max(3, round(span))}, matching \code{charLowess.m} lines 63-68.
#' @param iter  Number of bisquare robustness passes after the initial fit.
#'   \itemize{
#'     \item \code{iter = 0}: plain lowess (1 total WLS pass, no robustness).
#'     \item \code{iter = 4}: robust lowess, matching MATLAB's \code{'rlowess'}
#'       which uses \code{nIter = 5} total iterations (1 initial + 4 updates).
#'   }
#'
#' @return Numeric vector of smoothed values, same length as \code{y}.
#'
#' @details
#'   ## Algorithm
#'   For each point \eqn{i}, a window of exactly \eqn{k} points is selected:
#'   the window is centred on \eqn{i} at interior points and SHIFTED (not
#'   shrunk) at record boundaries, maintaining \eqn{k} neighbours throughout.
#'   Each neighbour \eqn{j} receives a tricubic weight:
#'   \deqn{w_j = \left(1 - \left(\frac{|i-j|}{d_{\max}}\right)^3\right)^3}
#'   where \eqn{d_{\max}} is the distance to the \strong{furthest} point in
#'   the window (not the half-window radius).  A weighted least-squares line
#'   is fitted and evaluated at \eqn{i}.  This matches
#'   \code{charLowess.m} lines 109-147 exactly.
#'
#'   For \code{iter > 0} (\code{'rlowess'}), bisquare robustness weights are
#'   computed from the residuals after each WLS pass and multiplied into the
#'   tricubic weights for the next pass (matching \code{charLowess.m} lines
#'   151-157).
#'
#'   ## Why not stats::lowess()?
#'   \code{stats::lowess()} trims (shrinks) the window at boundaries and
#'   normalises tricubic weights by the half-window radius \eqn{hw}.  The two
#'   differences compound to ~1e-3 absolute error in \code{charBkg} for the CO
#'   dataset (500 yr / 15 yr = 33-point window, 17-point boundary zone).
#'
#' @seealso [char_smooth()], [char_thresh_local()]
char_lowess <- function(y, x = NULL, span = 0.1, iter = 0L) {

  n <- length(y)
  if (is.null(x)) x <- seq_len(n)

  # ---- Convert span to integer window width ---------------------------------
  # Mirrors charLowess.m lines 63-68:
  #   if span < 1:  k = max(3, round(span * n))   % fraction -> point count
  #   else:         k = max(3, round(span))        % already a point count
  if (span < 1) {
    k <- max(3L, round(span * n))
  } else {
    k <- max(3L, round(span))
  }
  k  <- min(k, n)
  hw <- k %/% 2L                        # floor(k/2): nominal half-window

  n_iter <- 1L + as.integer(iter)       # total passes: 1 plain, 5 robust

  ys <- y
  rw <- rep(1.0, n)                     # robustness weights (all 1 initially)

  for (it in seq_len(n_iter)) {

    for (i in seq_len(n)) {

      # ---- Shifted window: always k points (charLowess.m lines 110-120) ----
      i0 <- i - hw
      i1 <- i0 + k - 1L
      if (i0 < 1L) { i0 <- 1L; i1 <- k        }
      if (i1 > n)  { i1 <- n;  i0 <- n - k + 1L }

      xi <- x[i0:i1]
      yi <- y[i0:i1]
      ri <- rw[i0:i1]

      # ---- Tricubic weights (charLowess.m lines 128-130) --------------------
      # Normalised by distance to the FURTHEST point in the window, not hw.
      d    <- abs(xi - x[i])
      dmax <- max(d) + .Machine$double.eps
      tri  <- (1.0 - (d / dmax)^3)^3
      w    <- tri * ri

      # ---- Weighted least-squares linear fit (charLowess.m lines 133-147) --
      sw   <- sum(w)
      swx  <- sum(w * xi)
      swy  <- sum(w * yi)
      swx2 <- sum(w * xi^2)
      swxy <- sum(w * xi * yi)
      det  <- sw * swx2 - swx^2

      if (abs(det) < .Machine$double.eps * max(abs(c(sw * swx2, swx^2)))) {
        ys[i] <- swy / (sw + .Machine$double.eps)   # near-singular: mean
      } else {
        a     <- (swy * swx2 - swx * swxy) / det
        b     <- (sw  * swxy - swx * swy)  / det
        ys[i] <- a + b * x[i]
      }
    }

    # ---- Bisquare robustness weights (charLowess.m lines 151-157) ----------
    if (iter > 0L && it < n_iter) {
      res <- y - ys
      s   <- stats::median(abs(res))
      if (s < .Machine$double.eps) break
      u   <- res / (6.0 * s)
      rw  <- pmax(0.0, 1.0 - u^2)^2 * (abs(u) < 1.0)
    }
  }

  ys
}
