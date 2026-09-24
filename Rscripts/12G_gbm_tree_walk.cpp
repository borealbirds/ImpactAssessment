// ---
// title: Impact Assessment: bit-identical replacement for gbm's compiled tree walk
// author: Mannfred Boehm
// ---
// Compiled by Rcpp::sourceCpp() from 12D; wrapped by 12G_gbm_tree_walk.R.
//
// gbm's gbm_pred (gbmentry.cpp) makes three non-inlined R API calls (REAL, ISNA, INTEGER)
// at EVERY node it visits. This does the same walk with the pointers hoisted, the NA test
// inlined and the trees flattened once per call: 3.6-4.4x faster on CAWA can11's models.
//
// The arithmetic is gbm_pred's exactly: out[i] = initF, then += the leaf of tree 0, 1, 2 ...
// in that order. Blocking only changes which observation is visited when, never the order of
// one observation's additions, so the sum is bit-identical (checked with identical() on
// 30-40k real pixels x 9150 trees, with and without missing values).

#include <Rcpp.h>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// R_IsNA, inlined: a NaN whose low word is 1954. gbm_pred sends ONLY this down the missing
// branch; any other NaN fails the `<` test and goes right. Replicated as-is.
static inline bool is_R_NA(double x) {
  if (!std::isnan(x)) return false;
  uint64_t u;
  std::memcpy(&u, &x, 8);
  return (uint32_t)(u & 0xFFFFFFFFu) == 1954u;
}

// [[Rcpp::export]]
NumericVector gbm_tree_walk(NumericMatrix X, List trees, List csplits, IntegerVector vartype,
                            double initF, int ntrees, int block) {
  const int n = X.nrow();
  if (ntrees > trees.size()) stop("ntrees exceeds the number of fitted trees");

  // flatten the trees: node arrays with global child offsets; kind -1 leaf, 0 cont, 1 cat
  std::vector<int> root(ntrees + 1);
  int total = 0;
  for (int t = 0; t < ntrees; t++) {
    root[t] = total;
    total += Rf_length(VECTOR_ELT(VECTOR_ELT(trees, t), 0));
  }
  std::vector<int> var(total), L(total), R(total), Mi(total), kind(total);
  std::vector<double> code(total);
  const int* vt = INTEGER(vartype);
  for (int t = 0; t < ntrees; t++) {
    SEXP tr = VECTOR_ELT(trees, t);
    const int* sv    = INTEGER(VECTOR_ELT(tr, 0));
    const double* sc = REAL(VECTOR_ELT(tr, 1));
    const int* le    = INTEGER(VECTOR_ELT(tr, 2));
    const int* ri    = INTEGER(VECTOR_ELT(tr, 3));
    const int* mi    = INTEGER(VECTOR_ELT(tr, 4));
    const int k = Rf_length(VECTOR_ELT(tr, 0)), r = root[t];
    for (int j = 0; j < k; j++) {
      var[r + j]  = sv[j];
      code[r + j] = sc[j];
      L[r + j] = r + le[j]; R[r + j] = r + ri[j]; Mi[r + j] = r + mi[j];
      kind[r + j] = (sv[j] == -1) ? -1 : (vt[sv[j]] == 0 ? 0 : 1);
    }
  }
  const int ncs = csplits.size();
  std::vector<int> csoff(ncs + 1, 0), csval;
  for (int c = 0; c < ncs; c++) {
    IntegerVector cv = as<IntegerVector>(csplits[c]);
    csoff[c] = csval.size();
    for (int j = 0; j < cv.size(); j++) csval.push_back(cv[j]);
  }
  csoff[ncs] = csval.size();

  const double* x = REAL(X);
  NumericVector out(n);
  double* o = REAL(out);
  for (int i = 0; i < n; i++) o[i] = initF;
  if (block <= 0) block = n;
  for (int b0 = 0; b0 < n; b0 += block) {
    const int b1 = std::min(n, b0 + block);
    for (int t = 0; t < ntrees; t++) {
      const int r = root[t];
      for (int i = b0; i < b1; i++) {
        int nd = r;
        while (kind[nd] != -1) {
          const double dx = x[(size_t)var[nd] * n + i];
          if (is_R_NA(dx)) nd = Mi[nd];
          else if (kind[nd] == 0) nd = (dx < code[nd]) ? L[nd] : R[nd];
          else {
            // gbm_pred indexes the split table with (int)dx unchecked, so a NaN or an
            // out-of-range level is undefined behaviour there (it segfaulted on a raw NaN
            // SurfaceWater_1km). Fail loudly instead; every in-range value is unchanged.
            const int c = (int)code[nd];
            if (std::isnan(dx) || dx < 0 || (int)dx >= csoff[c + 1] - csoff[c])
              stop("categorical value outside the model's levels (NaN or out of range)");
            const int ind = csval[csoff[c] + (int)dx];
            nd = (ind == -1) ? L[nd] : (ind == 1 ? R[nd] : Mi[nd]);
          }
        }
        o[i] += code[nd];
      }
    }
  }
  return out;
}

// rowsum(sc[kr, ], zones[kr], reorder = TRUE) placed on n_sub subbasin rows (a subbasin
// with no kept row stays 0), without materialising sc[kr, ]: 12F used to copy that once
// per coalition, 255 times per bootstrap, and the garbage helped OOM C1's workers.
// R's rowsum (src/main/unique.c) zeroes its result and adds x[j, col] in row order for
// each column; this does the same additions in the same order, so it is bit-identical.
// kr: 1-based, ascending rows of sc. zi: 1-based subbasin row of every row of sc, NA if none.
// [[Rcpp::export]]
NumericMatrix coal_rowsum(NumericMatrix sc, IntegerVector kr, IntegerVector zi, int n_sub) {
  const int n = sc.nrow(), p = sc.ncol(), m = kr.size();
  if (zi.size() != n) stop("coal_rowsum: zi must have one entry per row of sc");
  const int* k = INTEGER(kr);
  const int* z = INTEGER(zi);
  for (int a = 0; a < m; a++) {
    if (k[a] < 1 || k[a] > n || (a > 0 && k[a] <= k[a - 1]))
      stop("coal_rowsum: kr must be ascending row numbers of sc");
    if (z[k[a] - 1] != NA_INTEGER && (z[k[a] - 1] < 1 || z[k[a] - 1] > n_sub))
      stop("coal_rowsum: zi out of range");
  }
  NumericMatrix out(n_sub, p);
  const double* x = REAL(sc);
  double* o = REAL(out);
  for (int j = 0; j < p; j++) {
    const double* xj = x + (size_t)j * n;
    double* oj = o + (size_t)j * n_sub;
    for (int a = 0; a < m; a++) {
      const int r = k[a] - 1, g = z[r];
      if (g != NA_INTEGER) oj[g - 1] += xj[r];
    }
  }
  return out;
}
