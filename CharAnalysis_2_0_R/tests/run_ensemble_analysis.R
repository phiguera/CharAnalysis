# run_ensemble_analysis.R
# -----------------------------------------------------------------------
# Reference run -> load ensemble -> nearest-neighbor peak matching.
# Saving is handled by the calling script (char_run_ensemble.R) or
# manually after sourcing standalone.
#
# Can be sourced standalone from a fresh session:
#   setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")
#   source("CharAnalysis_2_0_R/tests/run_ensemble_analysis.R")
#
# When sourced from char_run_ensemble.R, `out` and `ensemble` are
# already in memory and both load steps are skipped automatically.
# -----------------------------------------------------------------------

library(CharAnalysis)

# When sourced from char_run_ensemble.R, tests_dir is already set and
# all paths are absolute -- skip setwd(). Only applies when run standalone.
if (!exists("tests_dir"))
  setwd("C:/Users/philip.higuera/OneDrive - The University of Montana/1_phiguera/1_working/CharAnalysis")

# ---- 1. Reference run --------------------------------------------------
# Skipped if `out` already exists in the workspace (e.g., sourced from
# char_run_ensemble.R, or already run this session).
if (!exists("out")) {
  message("Loading package...")
  devtools::load_all("CharAnalysis_2_0_R")
  if (!exists("params_file")) {
    stop("'out' not found and 'params_file' not set.\n",
         "  Set params_file <- 'path/to/MySite_charParams.csv' before sourcing,\n",
         "  or source char_run_ensemble.R instead (which sets params_file).")
  }
  message("Running reference CharAnalysis...")
  out <- CharAnalysis(params_file, plots = FALSE)
} else {
  message("Using existing 'out' from workspace.")
}


# ---- 2. Load ensemble --------------------------------------------------
# Skipped if `ensemble` already exists in the workspace.
if (!exists("ensemble")) {
  ensemble_rds <- if (exists("tests_dir")) {
    file.path(tests_dir, sprintf("%s_ensemble_results.rds", out$site))
  } else {
    sprintf("CharAnalysis_2_0_R/tests/%s_ensemble_results.rds", out$site)
  }
  if (!file.exists(ensemble_rds)) {
    stop("Ensemble RDS not found: ", ensemble_rds,
         "\nRun char_run_ensemble.R first to generate it.")
  }
  message("Loading ensemble from RDS...")
  ensemble <- readRDS(ensemble_rds)
} else {
  message("Using existing 'ensemble' from workspace.")
}

message(sprintf("  %d iterations  |  %d time steps  |  age grid %.0f-%.0f cal yr BP",
                ensemble$n_iter, length(ensemble$age_grid),
                min(ensemble$age_grid), max(ensemble$age_grid)))


# ---- 3. Reference peak ages --------------------------------------------
ccp           <- out$post$CharcoalCharPeaks
sel_col       <- ncol(ccp)
x             <- out$charcoal$ybpI
ref_peak_idx  <- which(ccp[, sel_col] > 0)
ref_peak_ages <- x[ref_peak_idx]
n_peaks       <- length(ref_peak_ages)



# ---- 4. Matching window (adaptive, per-peak) ---------------------------
# The window floor is mFRI_z / 2 for the zone containing each peak, where
# mFRI_z is the mean fire-return interval for that zone from the benchmark
# run, or yrInterp, whichever is larger. On top of that, the per-peak
# window expands where chronological uncertainty is high: ci95_width is
# interpolated from ensemble$chron_ci at the peak's age, and the window is
# max(mFRI_floor_z, ci95_width / 2). Using per-zone mFRI ensures the floor
# reflects local fire frequency rather than a record-wide average.
# Falls back to a 500-yr global default for zones with invalid mFRI, and
# to the zone floor for all peaks if chron_ci is unavailable.

yrInterp   <- out$pretreatment$yrInterp
zone_div   <- out$pretreatment$zoneDiv   # defined here; reused in section 9
n_zones    <- length(zone_div) - 1L
fri_params <- out$post$FRI_params_zone

# Per-zone mFRI floor: max(mFRI_z / 2, yrInterp).
mFRI_floor_by_zone <- vapply(seq_len(n_zones), function(z) {
  mz <- fri_params[z, 2L]
  if (!is.na(mz) && mz > 0) max(mz / 2, yrInterp) else 500
}, numeric(1L))

# Helper: return the mFRI floor for the zone containing a given age.
mFRI_floor_for_age <- function(age) {
  for (z in seq_len(n_zones)) {
    zb <- sort(c(zone_div[z], zone_div[z + 1L]))
    if (age >= zb[1L] && age <= zb[2L]) return(mFRI_floor_by_zone[z])
  }
  # Fallback: assign to nearest zone boundary
  mFRI_floor_by_zone[which.min(
    sapply(seq_len(n_zones), function(z) {
      zb <- sort(c(zone_div[z], zone_div[z + 1L]))
      min(abs(age - zb))
    })
  )]
}

# Per-peak mFRI floor (one value per reference peak).
mFRI_floor_k <- vapply(ref_peak_ages, mFRI_floor_for_age, numeric(1L))

# Interpolate ci95_width at each reference peak age from chron_ci.
chron_ci_ok <- exists("ensemble") && !is.null(ensemble$chron_ci)
if (chron_ci_ok) {
  ci95_at_peak    <- approx(
    x    = ensemble$chron_ci$median_age,
    y    = ensemble$chron_ci$ci95_width,
    xout = ref_peak_ages,
    rule = 2L
  )$y
  match_halfwin_k <- pmax(mFRI_floor_k, ci95_at_peak / 2)
} else {
  match_halfwin_k <- mFRI_floor_k
}

# Keep a scalar for orphan-stage checks (use median window).
match_halfwin <- median(match_halfwin_k)


# ---- 5. Nearest-neighbor matching --------------------------------------
age_grid  <- ensemble$age_grid
peaks_mat <- ensemble$peaks_matrix
n_iter    <- ensemble$n_iter

assigned       <- vector("list", n_peaks)
for (k in seq_len(n_peaks)) assigned[[k]] <- numeric(0)

# matched_matrix tracks which grid cells were claimed by a reference peak.
# Cells left at 0 after the loop are "orphan" detections -- ensemble peaks
# with no corresponding reference peak.
matched_matrix <- matrix(0L, nrow = nrow(peaks_mat), ncol = ncol(peaks_mat))

for (k in seq_len(n_peaks)) {
  win_k <- match_halfwin_k[k]
  for (i in seq_len(n_iter)) {
    det_idx <- which(peaks_mat[, i] == 1L)
    if (length(det_idx) == 0L) next
    det_ages <- age_grid[det_idx]
    dists    <- abs(det_ages - ref_peak_ages[k])
    best     <- which.min(dists)
    if (dists[best] <= win_k) {
      assigned[[k]] <- c(assigned[[k]], det_ages[best])
      matched_matrix[det_idx[best], i] <- 1L   # mark this cell as claimed
    }
  }
}


# ---- 6. Per-peak summaries (reference peaks) ---------------------------
make_peak_summary <- function(assigned_list, ref_ages, n_iter_val) {
  data.frame(
    ref_age    = ref_ages,
    n_iter     = n_iter_val,
    prob       = sapply(assigned_list, function(v) length(v) / n_iter_val),
    n_detected = sapply(assigned_list, length),
    mean_age   = sapply(assigned_list, function(v) if (length(v) > 0L) mean(v)            else NA_real_),
    median_age = sapply(assigned_list, function(v) if (length(v) > 0L) median(v)          else NA_real_),
    sd_age     = sapply(assigned_list, function(v) if (length(v) > 1L) sd(v)              else NA_real_),
    ci80_lo    = sapply(assigned_list, function(v) if (length(v) > 0L) quantile(v, 0.10)  else NA_real_),
    ci80_hi    = sapply(assigned_list, function(v) if (length(v) > 0L) quantile(v, 0.90)  else NA_real_),
    ci95_lo    = sapply(assigned_list, function(v) if (length(v) > 0L) quantile(v, 0.025) else NA_real_),
    ci95_hi    = sapply(assigned_list, function(v) if (length(v) > 0L) quantile(v, 0.975) else NA_real_)
  )
}

peak_summaries               <- make_peak_summary(assigned, ref_peak_ages, n_iter)
peak_summaries$type          <- "reference"
peak_summaries$proximity     <- NA_character_   # not applicable for reference peaks
peak_summaries$match_halfwin <- round(match_halfwin_k, 1)

any_over_1 <- sum(peak_summaries$prob > 1)
if (any_over_1 > 0L) {
  warning(sprintf("%d peaks have prob > 1.0 -- check matching window", any_over_1))
}



# ---- 7. Ensemble-only peak detection -----------------------------------
# Cells in peaks_mat that are 1 but were NOT claimed by any reference peak
# (matched_matrix == 0) are "orphan" detections. Contiguous clusters of
# these on the common_grid above `orphan_thresh` are candidate events that
# the single best-estimate chronology did not detect.
#
# Terminology: "ensemble-only peaks" -- neutral; avoids asserting that the
# reference run produced false negatives, since that requires additional
# evidence beyond detection frequency alone.
#
# Each orphan cluster uses the same nearest-neighbor assignment as
# reference peaks: for each cluster and each iteration, find the closest
# unmatched detection within match_halfwin and record its age. This gives
# per-event CI and detection frequency comparable to reference peaks.
# Note: two orphan clusters within 2 * match_halfwin of each other could
# share detections; this is flagged but not further resolved here.

orphan_thresh  <- 0.10    # min fraction of iterations to call a cluster a candidate

orphan_matrix  <- peaks_mat
orphan_matrix[matched_matrix == 1L] <- 0L
prob_orphan    <- rowMeans(orphan_matrix, na.rm = TRUE)

is_cand <- prob_orphan >= orphan_thresh

if (any(is_cand, na.rm = TRUE)) {

  # Find contiguous runs above threshold
  rl      <- rle(is_cand)
  ends    <- cumsum(rl$lengths)
  starts  <- ends - rl$lengths + 1L

  orphan_rep_ages  <- numeric(0)
  orphan_assgn     <- list()
  orphan_halfwins  <- numeric(0)

  for (r in which(rl$values)) {
    seg_idx  <- starts[r]:ends[r]
    peak_pos <- seg_idx[which.max(prob_orphan[seg_idx])]
    peak_age <- age_grid[peak_pos]

    # Adaptive window for this orphan cluster (same logic as reference peaks,
    # using the zone-specific mFRI floor for the age of this cluster).
    win_orp <- if (chron_ci_ok) {
      ci95_here <- approx(
        x    = ensemble$chron_ci$median_age,
        y    = ensemble$chron_ci$ci95_width,
        xout = peak_age,
        rule = 2L
      )$y
      max(mFRI_floor_for_age(peak_age), ci95_here / 2)
    } else {
      mFRI_floor_for_age(peak_age)
    }

    # Collect unmatched detected ages within win_orp for this cluster
    ages_in_win <- numeric(0)
    for (i in seq_len(n_iter)) {
      unm_idx <- which(orphan_matrix[, i] == 1L)
      if (length(unm_idx) == 0L) next
      unm_ages <- age_grid[unm_idx]
      dists    <- abs(unm_ages - peak_age)
      best     <- which.min(dists)
      if (dists[best] <= win_orp) {
        ages_in_win <- c(ages_in_win, unm_ages[best])
      }
    }

    orphan_rep_ages <- c(orphan_rep_ages, peak_age)
    orphan_assgn    <- c(orphan_assgn, list(ages_in_win))
    orphan_halfwins <- c(orphan_halfwins, win_orp)
  }

  orphan_summaries               <- make_peak_summary(orphan_assgn, orphan_rep_ages, n_iter)
  orphan_summaries$type          <- "ensemble_only"
  orphan_summaries$match_halfwin <- round(orphan_halfwins, 1)

  # Proximity flag: does this orphan fall within the adaptive window of any
  # reference peak? "near_reference" = yes (secondary detection near a known
  # event); "independent" = no (no nearby reference peak within its window).
  # Near-reference orphans arise because the matching loop claims at most one
  # detection per reference peak per iteration; a second detection in the same
  # iteration near the same reference peak is left unclaimed even if it falls
  # within the window. Independent orphans have no corresponding reference peak.
  orphan_summaries$proximity <- sapply(orphan_rep_ages, function(oa) {
    dists <- abs(ref_peak_ages - oa)
    nearest <- which.min(dists)
    if (dists[nearest] <= match_halfwin_k[nearest]) "near_reference" else "independent"
  })

  # Keep only independent orphans; near-reference cases handled by Section 7b.
  orphan_summaries <- orphan_summaries[orphan_summaries$proximity == "independent", ]
  if (nrow(orphan_summaries) == 0L) orphan_summaries <- NULL
  n_indp    <- if (is.null(orphan_summaries)) 0L else nrow(orphan_summaries)
  n_orphans <- n_indp

  # Warn if any two independent orphan clusters have overlapping windows
  if (n_indp > 1L) {
    indp_ages <- orphan_summaries$ref_age
    indp_wins <- orphan_summaries$match_halfwin
    for (j1 in seq_len(n_indp - 1L)) {
      for (j2 in (j1 + 1L):n_indp) {
        if (abs(indp_ages[j1] - indp_ages[j2]) < indp_wins[j1] + indp_wins[j2]) {
          message(sprintf(
            "  NOTE: independent orphan clusters at %.0f and %.0f cal yr BP are within\n        each other's search windows and may represent a single shifted event.",
            indp_ages[j1], indp_ages[j2]))
        }
      }
    }
  }

} else {
  orphan_summaries <- NULL
  n_orphans        <- 0L
  n_indp           <- 0L
}


# ---- 7b. Secondary detections within reference windows -----------------
# Many ensemble iterations detect more peaks than the reference run. The
# matching loop claims only the NEAREST detection per reference peak per
# iteration; additional detections within the same window are left unclaimed
# in orphan_matrix. Because these secondary detections are scattered across
# different time steps by chronological uncertainty, they rarely cluster
# enough at any single time step to exceed orphan_thresh in the time-step-
# based orphan scan above.
#
# This section aggregates unclaimed detections at the WINDOW level: for
# each reference peak k and each iteration i, find unclaimed detections
# within +/- match_halfwin_k of ref_peak_ages[k] and take the one furthest
# from the reference peak age. If >= orphan_thresh of iterations contribute
# a secondary detection near reference peak k, it is flagged as an
# ensemble-only peak of proximity "near_reference".
#
# ref_age in the output is the median detected age of the secondary
# detections (not the reference peak age), showing where the secondary
# event actually clusters.

secondary_assgn    <- list()
secondary_rep_ages <- numeric(0)   # median detected age (reported as ref_age)
secondary_halfwins <- numeric(0)

for (k in seq_len(n_peaks)) {
  win_k  <- match_halfwin_k[k]
  age_lo <- ref_peak_ages[k] - win_k
  age_hi <- ref_peak_ages[k] + win_k
  in_win <- which(age_grid >= age_lo & age_grid <= age_hi)

  sec_ages <- numeric(0)
  for (i in seq_len(n_iter)) {
    unm_idx <- in_win[which(orphan_matrix[in_win, i] == 1L)]
    if (length(unm_idx) == 0L) next
    ages_unm <- age_grid[unm_idx]
    furthest <- ages_unm[which.max(abs(ages_unm - ref_peak_ages[k]))]
    sec_ages <- c(sec_ages, furthest)
  }

  if (length(sec_ages) / n_iter >= orphan_thresh) {
    secondary_assgn    <- c(secondary_assgn, list(sec_ages))
    secondary_rep_ages <- c(secondary_rep_ages, median(sec_ages))
    secondary_halfwins <- c(secondary_halfwins, win_k)
  }
}

# ---- 7b (merging pass): collapse overlapping secondary candidates -------
# Candidates from adjacent reference windows can overlap, producing multiple
# entries for what is ecologically a single secondary event. Sort by median
# age, then greedily merge any candidate whose median falls within mFRI_floor
# of the running group median. After merging, recompute the group median from
# all pooled detected ages.
# Note: if the same iteration contributed a detection to two overlapping
# windows, it will appear twice in the merged pool. Detection frequency is
# therefore a slight overestimate for merged events; it is capped at 100%.

if (length(secondary_assgn) > 0L) {

  # Sort oldest-to-youngest (decreasing age BP = older first)
  ord                <- order(secondary_rep_ages, decreasing = TRUE)
  secondary_assgn    <- secondary_assgn[ord]
  secondary_rep_ages <- secondary_rep_ages[ord]
  secondary_halfwins <- secondary_halfwins[ord]

  merged_assgn    <- list()
  merged_rep_ages <- numeric(0)
  merged_halfwins <- numeric(0)

  k <- 1L
  while (k <= length(secondary_assgn)) {
    grp_ages <- secondary_assgn[[k]]
    grp_hw   <- secondary_halfwins[k]

    # Absorb the next candidate if its median falls within the zone mFRI floor of
    # the current group's running median
    while (k + 1L <= length(secondary_assgn)) {
      if (abs(secondary_rep_ages[k + 1L] - median(grp_ages)) <= mFRI_floor_for_age(median(grp_ages))) {
        k        <- k + 1L
        grp_ages <- c(grp_ages, secondary_assgn[[k]])
        grp_hw   <- max(grp_hw, secondary_halfwins[k])
      } else {
        break
      }
    }

    merged_assgn    <- c(merged_assgn, list(grp_ages))
    merged_rep_ages <- c(merged_rep_ages, median(grp_ages))
    merged_halfwins <- c(merged_halfwins, grp_hw)
    k <- k + 1L
  }

  # Re-apply threshold after merging (merged pool can only grow, but check
  # ensures consistency if any group somehow fell below after deduplication)
  keep <- sapply(merged_assgn,
                 function(v) min(length(v), n_iter) / n_iter >= orphan_thresh)
  secondary_assgn    <- merged_assgn[keep]
  secondary_rep_ages <- merged_rep_ages[keep]
  secondary_halfwins <- merged_halfwins[keep]
}

if (length(secondary_assgn) > 0L) {
  secondary_summaries               <- make_peak_summary(secondary_assgn,
                                                         secondary_rep_ages, n_iter)
  # Cap prob at 1.0 in case of double-counting from overlapping windows
  secondary_summaries$prob          <- pmin(secondary_summaries$prob, 1.0)
  secondary_summaries$type          <- "ensemble_only"
  secondary_summaries$proximity     <- "near_reference"
  secondary_summaries$match_halfwin <- round(secondary_halfwins, 1)

  # Coherence filter: retain only candidates whose 95% CI spans < 50% of the
  # full matching window. ci_frac = ci95_width / (2 * match_halfwin).
  # Values near 1.0 indicate detections scattered across the entire window
  # (noise); values < 0.5 indicate genuine temporal clustering at a
  # consistent age (real secondary event).
  ci_frac <- (secondary_summaries$ci95_hi - secondary_summaries$ci95_lo) /
               (secondary_summaries$match_halfwin * 2)
  secondary_summaries <- secondary_summaries[ci_frac < 0.5, ]
  if (nrow(secondary_summaries) == 0L) secondary_summaries <- NULL

  n_secondary <- if (is.null(secondary_summaries)) 0L else nrow(secondary_summaries)
} else {
  secondary_summaries <- NULL
  n_secondary         <- 0L
}


# ---- 8. Build output data frame (charResults-style column names) -------
# Combines reference peaks and ensemble-only peaks into a single CSV.
# type column: "reference" | "ensemble_only"
# Kept in memory as `out_df`; saving handled by char_run_ensemble.R.

build_out_df <- function(ps) {
  data.frame(
    "age Top_i (yr BP)"       = ps$ref_age,
    "type"                    = ps$type,
    "proximity"               = if ("proximity" %in% names(ps)) ps$proximity else NA_character_,
    "match halfwin (yr)"      = ps$match_halfwin,
    "n ensemble (iterations)" = ps$n_iter,
    "detection freq (%)"      = round(ps$prob * 100, 1),
    "n detected (peaks)"      = ps$n_detected,
    "mean age (yr BP)"        = round(ps$mean_age,   0),
    "median age (yr BP)"      = round(ps$median_age, 0),
    "sd age (yr)"             = round(ps$sd_age,     1),
    "CI80 lo (yr BP)"         = round(ps$ci80_lo,    1),
    "CI80 hi (yr BP)"         = round(ps$ci80_hi,    1),
    "CI95 lo (yr BP)"         = round(ps$ci95_lo,    1),
    "CI95 hi (yr BP)"         = round(ps$ci95_hi,    1),
    check.names = FALSE
  )
}

all_summaries <- rbind(peak_summaries, orphan_summaries, secondary_summaries)
# Sort young-to-old (smallest age BP at top, consistent with charResults convention)
all_summaries <- all_summaries[order(all_summaries$ref_age, decreasing = FALSE), ]

out_df <- build_out_df(all_summaries)



# ---- 9. Summary statistics ----------------------------------------------
# Two quantities reported separately (not combined):
#
#   (a) Peak timing uncertainty -- SD of detected age across ensemble
#       iterations. Reflects how much chronological uncertainty shifts the
#       absolute timing of individual fire events. Units: years.
#
#   (b) Mean FRI uncertainty -- distribution of arithmetic mean FRI values
#       across ensemble iterations, for the whole record and by zone.
#       Reported as: min | median | max across iterations, alongside the
#       reference-run arithmetic mean FRI for direct comparison.
#       Arithmetic mean used throughout (not Weibull mean) because the
#       Weibull fit is unreliable for short records or zones with few peaks.
#
# NOT combined with ecological variability in FRIs. Ecological variability
# (the irregular spacing of actual fire events) is captured by the Weibull
# shape parameter (k) from the reference run and is typically much larger
# than chronological uncertainty. Combining the two would misrepresent the
# fire-regime signal by inflating apparent timing uncertainty.

# Helper: arithmetic mean FRI from a vector of peak ages (cal yr BP).
# Sorts descending (oldest first), then computes FRI for each consecutive
# pair as older_age - younger_age (= positive interval in years), starting
# from the second oldest peak. This is -diff() of the descending sequence.
arith_mean_fri <- function(ages) {
  ages_sorted <- sort(ages, decreasing = TRUE)   # oldest (largest cal yr BP) first
  if (length(ages_sorted) < 2L) return(NA_real_)
  mean(-diff(ages_sorted))                        # older - younger = positive FRI
}

message("\n", paste(rep("=", 60), collapse = ""))
message("  CHRONOLOGICAL UNCERTAINTY SUMMARY")
message(paste(rep("=", 60), collapse = ""))

# Ensemble overview
message(sprintf(
  "  Site         : %s
  Ensemble     : %d iterations  |  %.0f - %.0f cal yr BP
  Reference run: %d peaks  |  %.0f - %.0f cal yr BP",
  out$site,
  ensemble$n_iter,
  min(ensemble$age_grid), max(ensemble$age_grid),
  n_peaks,
  min(ref_peak_ages), max(ref_peak_ages)
))

# -- (a) Detection frequency and matching window -------------------------
det_freq_pct <- peak_summaries$prob * 100
hw           <- peak_summaries$match_halfwin   # half-window (± yr) per peak
message(sprintf(
  "\n(a) Detection frequency (%% of %d iterations):\n      Median: %.1f%%  |  range: %.1f - %.1f%%  (N = %d peaks)\n      >= 90%%: %d peaks  |  >= 75%%: %d peaks  |  >= 50%%: %d peaks\n    Matching window (+/- yr around each reference peak age):\n      Median: +/- %.0f yr  |  range: +/- %.0f to +/- %.0f yr\n      Note: window = max(mFRI_z / 2, yrInterp, ci95_width / 2) per zone; wide windows\n            may absorb potential ensemble-only peaks near reference peaks.",
  ensemble$n_iter,
  median(det_freq_pct), min(det_freq_pct), max(det_freq_pct), n_peaks,
  sum(det_freq_pct >= 90), sum(det_freq_pct >= 75), sum(det_freq_pct >= 50),
  median(hw), min(hw), max(hw)
))

# -- (b) Peak timing SD --------------------------------------------------
ref_sd <- peak_summaries$sd_age[!is.na(peak_summaries$sd_age)]
message(sprintf(
  "\n(b) Reference peak timing uncertainty (1 SD, chronological uncertainty only)
      Median SD : %.0f yr  |  range: %.0f - %.0f yr  (N = %d peaks)",
  median(ref_sd), min(ref_sd), max(ref_sd), length(ref_sd)
))

# Per-zone SD summary if multiple zones exist
if (n_zones > 1L) {
  message("      By zone:")
  for (z in seq_len(n_zones)) {
    z_bounds  <- sort(c(zone_div[z], zone_div[z + 1L]))
    in_zone_k <- which(ref_peak_ages >= z_bounds[1L] & ref_peak_ages <= z_bounds[2L])
    sd_z      <- peak_summaries$sd_age[in_zone_k]
    sd_z      <- sd_z[!is.na(sd_z)]
    if (length(sd_z) == 0L) {
      message(sprintf("        Zone %d: no peaks with SD estimate", z))
    } else {
      message(sprintf(
        "        Zone %d: median %.0f yr | range %.0f - %.0f yr | n = %d",
        z, median(sd_z), min(sd_z), max(sd_z), length(sd_z)))
    }
  }
}

# -- 9b. Ensemble mean FRI -----------------------------------------------
# For each ensemble iteration, compute the arithmetic mean inter-event
# interval from peaks_matrix using arith_mean_fri() (descending sort,
# -diff()). Requires >= 2 detected peaks per zone/record to contribute
# a FRI estimate.

# Peaks per iteration summary
peaks_per_iter <- colSums(peaks_mat, na.rm = TRUE)
message(sprintf(
  "\n(c) Peaks detected per ensemble iteration:\n      Reference run: %d peaks\n      Ensemble -- median: %.0f  |  range: %.0f - %.0f  (N = %d iterations)\n      Note: ref and ensemble mean FRI are not directly comparable if peak\n            counts differ -- chronological uncertainty affects detection,\n            not just timing.",
  n_peaks,
  median(peaks_per_iter), min(peaks_per_iter), max(peaks_per_iter),
  n_iter
))

message("\n(d) Mean FRI (chronological uncertainty | reference-run arithmetic mean)")
message("    Ensemble range reflects chronological uncertainty only.")

fri_ensemble <- vector("list", n_zones + 1L)  # index 1 = whole record, 2:(n_zones+1) = zones
names(fri_ensemble) <- c("whole_record", paste0("zone_", seq_len(n_zones)))

# Whole-record reference arithmetic mean FRI
ref_mean_fri_overall <- arith_mean_fri(ref_peak_ages)

# Whole-record ensemble FRI
fri_iter_overall <- vapply(seq_len(n_iter), function(i) {
  arith_mean_fri(age_grid[peaks_mat[, i] == 1L])
}, numeric(1L))
fri_iter_overall_valid  <- fri_iter_overall[!is.na(fri_iter_overall)]
fri_ensemble[["whole_record"]] <- fri_iter_overall_valid

message(sprintf(
  "    Whole record | ref mean FRI: %.0f yr | ensemble: min %.0f  median %.0f  max %.0f yr  (n_iter = %d)",
  ref_mean_fri_overall,
  min(fri_iter_overall_valid), median(fri_iter_overall_valid),
  max(fri_iter_overall_valid), length(fri_iter_overall_valid)
))

# Zone-level FRI
if (n_zones > 1L) {
  for (z in seq_len(n_zones)) {
    z_bounds <- sort(c(zone_div[z], zone_div[z + 1L]))

    # Reference run arithmetic mean FRI for this zone
    in_zone_ref  <- ref_peak_ages >= z_bounds[1L] & ref_peak_ages <= z_bounds[2L]
    ref_mean_fri_z <- arith_mean_fri(ref_peak_ages[in_zone_ref])

    # Ensemble FRI for this zone
    fri_iter_z <- vapply(seq_len(n_iter), function(i) {
      arith_mean_fri(
        age_grid[peaks_mat[, i] == 1L &
                 age_grid >= z_bounds[1L] &
                 age_grid <= z_bounds[2L]]
      )
    }, numeric(1L))

    fri_iter_z_valid <- fri_iter_z[!is.na(fri_iter_z)]
    fri_ensemble[[paste0("zone_", z)]] <- fri_iter_z_valid
    n_valid <- length(fri_iter_z_valid)

    if (n_valid < 10L) {
      message(sprintf(
        "    Zone %d: fewer than 10 iterations with >= 2 peaks; FRI estimate unreliable.", z))
    } else {
      ref_str <- if (is.na(ref_mean_fri_z)) "NA" else sprintf("%.0f yr", ref_mean_fri_z)
      message(sprintf(
        "    Zone %d       | ref mean FRI: %s | ensemble: min %.0f  median %.0f  max %.0f yr  (n_iter = %d)",
        z, ref_str,
        min(fri_iter_z_valid), median(fri_iter_z_valid),
        max(fri_iter_z_valid), n_valid))
    }
  }
}

# -- 9c. Ensemble-only peaks ---------------------------------------------
n_ensemble_only <- n_secondary + n_indp
message(sprintf(
  "\n(e) Ensemble-only peaks (>= %.0f%% detection threshold): %d total",
  orphan_thresh * 100, n_ensemble_only
))

# Near-reference: window-level secondary detections within reference windows
if (n_secondary > 0L) {
  sec_ages_str <- paste(
    round(sort(secondary_summaries$ref_age, decreasing = TRUE)),
    collapse = ", "
  )
  message(sprintf(
    "    Near-reference: %d  |  median ages (cal yr BP): %s\n      Unclaimed detections within a reference peak's window, present in\n      >= %.0f%% of iterations, with 95%% CI spanning < 50%% of window (ci_frac < 0.5).",
    n_secondary, sec_ages_str, orphan_thresh * 100))
} else {
  message(sprintf(
    "    Near-reference: 0  (no reference window had a secondary cluster\n      in >= %.0f%% of iterations)",
    orphan_thresh * 100))
}

# Independent: time-step orphan clusters outside all reference windows
if (n_indp > 0L) {
  indp_ages_out <- sort(orphan_summaries$ref_age, decreasing = TRUE)
  message(sprintf(
    "    Independent: %d  |  ages (cal yr BP): %s\n      Detected outside all reference windows in >= %.0f%% of iterations.",
    n_indp,
    paste(round(indp_ages_out), collapse = ", "),
    orphan_thresh * 100))
} else {
  message(sprintf(
    "    Independent: 0  (no clusters outside reference windows exceeded\n      %.0f%% detection threshold)",
    orphan_thresh * 100))
}

message(paste(rep("=", 60), collapse = ""))

# Store summary as a list for downstream access
ensemble_fri_summary <- list(
  peak_sd_median      = median(ref_sd),
  peak_sd_range       = range(ref_sd),
  ref_mean_fri        = ref_mean_fri_overall,
  fri_by_zone         = fri_ensemble
)
