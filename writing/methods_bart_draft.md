## Backfilling: concept

We used Bayesian Additive Regression Trees (BART; Chipman et al. 2010) to model the
relationship between biotic (vegetation) features and the abiotic environment, fitting a
separate model within each watershed using only that watershed's low-footprint pixels. BART
is a sum-of-trees model: rather than fitting one large regression tree, it fits an ensemble of
many deliberately shallow trees and adds their predictions together. Each tree is a weak
learner that captures only a small part of the signal, and the ensemble as a whole
accommodates non-linear responses and interactions among predictors without the analyst
specifying either in advance.

For a continuous biotic feature y at pixel i we modelled

    log(y_i + 1) = SUM_{j=1..m} g(x_i ; T_j , M_j) + e_i ,    e_i ~ Normal(0, sigma^2)      (1)

where m is the number of trees, T_j is the structure of the j-th tree (its binary splitting
rules), M_j is the set of values at that tree's terminal nodes, g(x_i ; T_j , M_j) returns the
terminal-node value of tree j for the pixel whose predictors are x_i, and sigma^2 is the
residual variance. The predictor vector x_i holds the abiotic covariates, the geodetic
coordinates, and the biotic features preceding y in the backfilling hierarchy. We modelled
log(y_i + 1) rather than y_i because cover variables are strongly right-skewed and bounded
below at zero.

There is no closed-form posterior for a model of this form, so it is fit by Markov chain Monte
Carlo: Bayesian backfitting (Hastie and Tibshirani 2000) implemented as a systematic-scan
Gibbs sampler (Geman and Geman 1984) in which trees are updated one at a time. At iteration h
the sampler visits tree T_j and computes its partial residuals — the response minus the summed
predictions of the other m - 1 trees, taking trees T_1 ... T_(j-1) as already updated within
iteration h and trees T_(j+1) ... T_m as they stood at iteration h - 1. Tree j is then fit to
those partial residuals alone, so each tree only ever explains what the rest of the ensemble
has left over.

Updating a tree takes two steps. First, a small structural change is proposed — growing a pair
of terminal nodes, pruning a pair back, or changing a splitting rule (Chipman et al. 2010) —
producing a candidate structure T_j*. The candidate is accepted or rejected by a
Metropolis-Hastings step whose acceptance probability is the ratio of the marginal likelihoods
of T_j* and T_j, multiplied by the ratio of their prior probabilities and by the ratio of the
proposal probabilities. Those marginal likelihoods are available in closed form because the
terminal-node values are integrated out analytically, which is what makes the sampler
tractable for ensembles of this size. Second, conditional on the accepted structure, new
terminal-node values M_j are drawn from their conjugate posterior. Once every tree has been
visited, sigma^2 is drawn from its own conjugate posterior, completing one iteration.
Repeating this scheme many times, and discarding an initial burn-in period, yields draws from
the joint posterior of the trees and sigma^2.

What this returns is not a single fitted relationship but a posterior distribution over
sum-of-trees functions, each of which can be evaluated at a pixel that was not used in
training. Prediction therefore yields a distribution of plausible values at each pixel rather
than a point estimate, and we carry that distribution — not only its mean — forward into the
bird density models (see XYZ), so that uncertainty in the backfilled vegetation propagates
into the counterfactual population estimates.

Because BART is non-parametric, the prior encodes no ecological hypothesis about how y
responds to x: no functional form, no specified interactions. What it does encode is a
preference for simplicity, through three components. (i) A tree-structure prior makes deep
trees unlikely: a node at depth d splits further with probability alpha(1 + d)^(-beta), with
d = 0 at the root, so splitting becomes rapidly less probable with depth. (ii) A prior on the
terminal-node values shrinks each tree's contribution toward zero, scaled so that the m trees
together span the observed range of the response; a shrinkage parameter k controls its width,
with larger k tightening the prior so that each individual tree contributes less. (iii) A
prior on sigma^2 anchored to a rough data-based estimate of the residual standard deviation,
which stops the ensemble from driving residual variance to zero by overfitting. Together these
make each tree a weak learner, so the fitted surface is assembled from many small
contributions rather than a few deep, confidently specified interactions. The implied
ecological assumption is correspondingly modest: that the biotic response can be built up from
many simple, shallow rules applied to the abiotic environment, with the data left to determine
which abiotic predictors matter.

## Backfilling: training

We implemented BART with the `BART` package v2.9.9 (Sparapani et al. 2021) in R v4.5.1 (R Core
Team 2025). We integrated Claude Opus 5 (Anthropic PBC, 2026) with RStudio v2025.9.1.401
(Posit Team, 2025) to co-develop the BART pipeline and subsequent counterfactual analyses. All
code was reviewed and edited by MMAB.

Models were fit independently within each watershed. Within a watershed, the training set was
its low-footprint pixels and the prediction targets were its high-footprint pixels. Latitude
and longitude were centred and scaled within each watershed and included as predictors, which
places them on a scale comparable to the other predictors and lets the trees absorb residual
spatial structure not explained by the abiotic covariates. BART admits no missing predictor
values — a single missing entry returns a missing prediction for that pixel — so where a
predictor was partially missing we replaced the missing entries with the training-set median
and added a binary missingness indicator as an additional predictor, allowing the model to
split on "was missing" where missingness is itself informative (e.g. phenology metrics over
open water). Random-number seeds were set deterministically from the watershed and covariate
indices, so every fit is reproducible.

Continuous biotic features were fit with `gbart()` (`type = "wbart"`) using m = 50 trees and
700 retained posterior draws following 300 burn-in iterations, with no thinning. We set the
terminal-node shrinkage parameter to k = 3 (package default 2), shrinking each tree's
contribution more strongly than the default — a deliberately conservative choice given that
these models are applied to pixels outside the training set. We retained the package defaults
for the tree-structure prior, alpha = 0.95 and beta = 2, so a node at depth d splits with
probability 0.95(1 + d)^(-2). The prior on the residual standard deviation was centred on the
standard deviation of the log-transformed training response, with 3 prior degrees of freedom.
We enabled the sparsity-inducing Dirichlet prior on splitting-variable probabilities
(`sparse = TRUE`; Linero 2018): rather than treating every predictor as equally eligible at
every split, the sampler learns a probability vector over predictors and concentrates splits
on the informative ones. This matters here because the abiotic predictor set is large relative
to the number of low-footprint training pixels in many watersheds. Responses were
back-transformed as exp(.) - 1 for reporting and downstream use. At every backfilled pixel we
retained 100 draws sampled at random from the 700-draw posterior and wrote them to disk; these
are resampled when re-predicting bird density (XYZ), so that BART's predictive uncertainty is
propagated rather than collapsed to a posterior mean.

Six categorical land-cover covariates (XYZ) were fit with `mbart2()` using binary probit
ensembles (`type = "pbart"`). This routine fits one ensemble per land-cover class, each
answering the binary question "is this pixel class h, or not". Within an ensemble the binary
response is handled by data augmentation: a continuous latent variable is drawn for each pixel
from a normal distribution truncated to be positive when the pixel belongs to the class and
negative when it does not (Albert and Chib 1993), so that the trees can then be updated
exactly as in the continuous case. The K ensemble scores at a pixel are exponentiated and
normalised to sum to one, giving a predicted probability for each class; we took the
highest-probability class as the backfilled land cover. Multinomial models used 40 trees,
k = 3, and 500 retained draws after 150 burn-in iterations, keeping every tenth draw. Thinning
is needed here but not for the continuous models because the latent-variable updates are
non-conjugate and consecutive draws are far more autocorrelated. Land-cover classes were
always backfilled last in the hierarchy, and categorical covariates were never used as
predictors of continuous biotic features.

We summarised within-sample performance per watershed x covariate as the root mean square
error (RMSE) and mean absolute error (MAE) of the posterior mean against observed values on
the original scale; the median and 95% credible interval of R^2, computed as 1 - SSE/SST
separately within each posterior draw; the posterior mean and standard deviation of the
residual standard deviation sigma (reported on the log(y + 1) scale on which it is estimated);
and 95% posterior predictive coverage, the proportion of observed values falling inside the
2.5-97.5% interval of the posterior predictive distribution. We assessed out-of-sample
performance by withholding a random 10% of each watershed's low-footprint pixels before
fitting and computing RMSE, MAE and R^2 on that holdout. For categorical covariates we report
classification accuracy and the mean Shannon entropy of the predicted class probabilities,
in-sample and on the holdout, together with confusion matrices (Tables XYZ).

A small number of watershed x covariate combinations do not support a fitted model: where a
biotic feature was constant across a watershed's low-footprint pixels, took two or fewer
distinct values after transformation, or where the watershed held a single high-footprint
pixel to backfill. In these cases we substituted the training-set mean (continuous features)
or modal class (categorical features) for a BART fit, and excluded the combination from the
performance tables.
