# ---
# title: Shapley value utilities for sector attribution
# author: Mannfred Boehm
# ---
# Enumerates all 2^N coalitions of N sectors, maps between coalition IDs
# and sector name vectors, and computes exact Shapley values from a table
# of coalition impact values.
#
# A coalition S is the set of sectors being "removed" (backfilled).
# v(S) = counterfactual population under S - observed population.
#
# Coalition IDs are 1-indexed integers in [1, 2^N].  ID 1 = empty coalition
# (no sectors removed, v = 0), ID 2^N = full coalition (all sectors removed).
# ---

#' Map coalition ID (1-indexed) to sector names
#' @param id Integer coalition ID (1 to 2^N)
#' @param sectors Character vector of all sector names (fixed order)
#' @return Character vector of sector names in this coalition
coalition_id_to_sectors <- function(id, sectors) {
  n <- length(sectors)
  bits <- as.logical(intToBits(id - 1L)[1:n])
  sectors[bits]
}

#' Map sector names to coalition ID (1-indexed)
#' @param coalition Character vector of sector names in this coalition
#' @param sectors Character vector of all sector names (fixed order)
#' @return Integer coalition ID
sectors_to_coalition_id <- function(coalition, sectors) {
  bits <- sectors %in% coalition
  sum(2L^(which(bits) - 1L)) + 1L
}

#' Total number of coalitions for N sectors
n_coalitions <- function(sectors) 2L^length(sectors)

#' Compute exact Shapley values from a named vector of coalition impacts
#'
#' @param v Named numeric vector: names are coalition IDs (as character strings),
#'   values are the population impact v(S) = cf(S) - obs for that coalition.
#'   Must include all 2^N coalitions (ID "1" for empty set should have v = 0).
#' @param sectors Character vector of all sector names (fixed, sorted order).
#' @return Named numeric vector of Shapley values (one per sector), summing
#'   to v(full coalition).
compute_shapley <- function(v, sectors) {
  n   <- length(sectors)
  phi <- setNames(numeric(n), sectors)

  for (j in seq_len(n)) {
    s_j    <- sectors[j]
    others <- sectors[-j]

    # iterate over all subsets S of N \ {j}
    for (bits in 0:(2^(n - 1) - 1)) {
      S      <- others[which(as.logical(intToBits(bits)[1:(n - 1)]))]
      s_size <- length(S)

      # Shapley weight
      w <- factorial(s_size) * factorial(n - s_size - 1L) / factorial(n)

      # marginal contribution: v(S u {j}) - v(S)
      id_with    <- sectors_to_coalition_id(c(S, s_j), sectors)
      id_without <- sectors_to_coalition_id(S, sectors)

      phi[j] <- phi[j] + w * (v[as.character(id_with)] - v[as.character(id_without)])
    }
  }

  phi
}

#' Shapley values as a linear map: phi = W %*% v
#'
#' Shapley values are linear in v, so compute_shapley() is a fixed [N x 2^N] matrix
#' applied to the vector of coalition impacts. Row j, column T (coalition ID) holds
#' +w(|T| - 1) if j is in T and -w(|T|) if not, with w(s) = s! (N - s - 1)! / N!.
#' Every column sums to 0 except the full coalition's, which sums to 1 (efficiency).
#' 12F applies it to each (bootstrap, scenario) sample, which a function of the mean
#' v cannot do.
#' @param sectors Character vector of all sector names (fixed, sorted order).
#' @return Numeric matrix [N sectors x 2^N coalition IDs], dimnames sectors x "1".."2^N".
shapley_weight_matrix <- function(sectors) {
  n <- length(sectors)
  w <- function(s) factorial(s) * factorial(n - s - 1L) / factorial(n)
  W <- matrix(0, n, 2L^n, dimnames = list(sectors, as.character(seq_len(2L^n))))
  for (id in seq_len(2L^n)) {
    in_T <- sectors %in% coalition_id_to_sectors(id, sectors)
    s    <- sum(in_T)
    if (s > 0L) W[in_T, id]  <- w(s - 1L)
    if (s < n)  W[!in_T, id] <- -w(s)
  }
  W
}

#' Canonical sector list used throughout the pipeline
#' Order matters: coalition IDs depend on this fixed ordering.
canonical_sectors <- function() {
  sort(c("built", "crop", "dam_and_associated_reservoir",
         "mines", "oil_gas", "pasture", "rail", "roads"))
}
