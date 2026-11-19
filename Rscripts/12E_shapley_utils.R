# ---
# title: Shapley value utilities for sector attribution
# author: Mannfred Boehm
# ---

# a coalition S is a set of sectors removed (backfilled); its impact v(S) = counterfactual
# population under S - observed population. coalition ids run 1..2^N for N sectors:
# id - 1 in binary says which sectors are in (id 1 = none, v = 0; id 2^N = all)

# the 8 sectors, in the fixed order the coalition ids depend on
canonical_sectors <- function() {
  sort(c("built", "crop", "dam_and_associated_reservoir", "mines", "oil_gas", "pasture", "rail", "roads"))
}

coalition_id_to_sectors <- function(id, sectors) {
  sectors[as.logical(intToBits(id - 1L)[seq_along(sectors)])]
}

sectors_to_coalition_id <- function(coalition, sectors) {
  sum(2L^(which(sectors %in% coalition) - 1L)) + 1L
}

# exact Shapley values from v, a vector of impacts named by coalition id ("1".."2^N"):
# phi_j = sum over S not containing j of |S|! (N - |S| - 1)! / N! * (v(S + j) - v(S)).
# they sum to v(all sectors)
compute_shapley <- function(v, sectors) {
  n   <- length(sectors)
  phi <- setNames(numeric(n), sectors)
  for (j in seq_len(n)) {
    others <- sectors[-j]
    for (bits in 0:(2^(n - 1) - 1)) {
      S <- others[which(as.logical(intToBits(bits)[1:(n - 1)]))]
      w <- factorial(length(S)) * factorial(n - length(S) - 1L) / factorial(n)
      phi[j] <- phi[j] + w * (v[as.character(sectors_to_coalition_id(c(S, sectors[j]), sectors))] -
                              v[as.character(sectors_to_coalition_id(S, sectors))])
    }
  }
  phi
}

# the same as a linear map, phi = W %*% v, with W [N sectors x 2^N coalition ids]: entry
# (j, T) is +w(|T| - 1) if j is in T and -w(|T|) if not, w(s) = s! (N - s - 1)! / N!. being
# linear, it applies to each (bootstrap, scenario) sample of v (12F), which compute_shapley()
# of the mean v cannot do
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
