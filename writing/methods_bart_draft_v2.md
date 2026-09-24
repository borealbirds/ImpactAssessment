## Backfilling: concept

Within each watershed we modelled each biotic feature as a function of the abiotic environment
using Bayesian Additive Regression Trees (BART; Chipman et al. 2010), fitted only to that
watershed's low-footprint pixels. BART represents the response as a sum of many shallow
regression trees. Each tree captures a small part of the signal, and together the ensemble
accommodates non-linear responses and interactions among predictors without the analyst
specifying either in advance.

For a continuous biotic feature y at pixel i,

    log(yᵢ + 1) = Σⱼ₌₁ᵐ g(xᵢ ; Tⱼ, Mⱼ) + εᵢ ,    εᵢ ~ N(0, σ²)        (1)

where m is the number of trees, Tⱼ is the structure of tree j (its binary splitting rules), Mⱼ
is the set of values at its terminal nodes, g(xᵢ ; Tⱼ, Mⱼ) returns the terminal-node value of
tree j for a pixel with predictors xᵢ, and σ² is the residual variance. The predictor vector xᵢ
holds the abiotic covariates, the pixel's projected coordinates, and the biotic features that
precede y in the backfilling hierarchy. We modelled log(y + 1) because vegetation features are
right-skewed and bounded below at zero.

The posterior of (1) was sampled by Markov chain Monte Carlo using Bayesian backfitting (Hastie
and Tibshirani 2000), a systematic-scan Gibbs sampler (Geman and Geman 1984) that updates the
trees one at a time. When the sampler reaches tree j in iteration h, it computes the partial
residuals: the response minus the summed predictions of the other m − 1 trees, taking trees 1
to j − 1 as already updated in iteration h and trees j + 1 to m as they stood at the end of
iteration h − 1. Tree j is then refitted to those partial residuals alone, so each tree explains
only what the rest of the ensemble has left over.

Each tree update has two steps. First, a local change to the tree's structure is proposed:
either a birth, which splits a terminal node into two, or a death, which collapses two sibling
terminal nodes back into their parent (Chipman et al. 2010; Sparapani et al. 2021). The proposed
structure Tⱼ* is accepted with Metropolis–Hastings probability min(1, r), where r is the product
of three ratios between Tⱼ* and Tⱼ: their marginal likelihoods, which measure how well each
structure explains the partial residuals; their prior probabilities, which penalise deeper
trees (see below); and their proposal probabilities, which correct for births and deaths not
being equally likely to be proposed. The marginal likelihoods have a closed form because the
normal prior on terminal-node values is conjugate to the normal likelihood, so node values can
be integrated out and a structure evaluated before any node values are chosen. Second,
conditional on the accepted structure, new terminal-node values Mⱼ are drawn from their normal
posterior. Once all m trees have been updated, σ² is drawn from its inverse-χ² posterior,
completing one iteration. After an initial burn-in is discarded, successive iterations are
draws from the joint posterior of the trees and σ².

The result is not a single fitted relationship but a posterior distribution over sum-of-trees
functions f(x) = Σⱼ g(x ; Tⱼ, Mⱼ), each of which can be evaluated at pixels outside the training
set. At a high-footprint pixel, the posterior draws of f(xᵢ) describe our uncertainty about the
expected value of the biotic feature given that pixel's abiotic environment. This is uncertainty
about the mean response, not a prediction of the value an individual pixel would take: it
excludes the residual variation (σ²) among pixels that share the same predictors, so a
backfilled surface is smoother than observed vegetation. We carried these draws, rather than
only their mean, forward into the bird density models (see XYZ), so that uncertainty in the
fitted vegetation–environment relationship propagates into the counterfactual population
estimates.

Because BART is non-parametric, its prior specifies no functional form and no interactions
between y and x. Instead it regularises model complexity through three components. (i) A
tree-structure prior under which a node at depth d (d = 0 at the root) splits with probability
α(1 + d)^−β, so splits become rapidly less likely as trees deepen. (ii) A normal prior on
terminal-node values, centred on zero and scaled so that the observed range of the training
response lies within ±k prior standard deviations of the sum of the m trees; larger k shrinks
each tree's contribution more strongly. (iii) A scaled inverse-χ² prior on σ², calibrated to a
rough data-based estimate of the residual standard deviation, which stops the ensemble from
overfitting by driving the residual variance towards zero. We also used a sparsity-inducing
prior on which predictors are chosen for splits (Linero 2018; see below). Together these make
each tree a weak learner. The implied ecological assumption is correspondingly modest: that a
biotic feature can be built up from many shallow, additive rules on the abiotic environment, and
that only some abiotic predictors are informative in a given watershed, with the data left to
determine which.

## Backfilling: training

We fitted BART with the `BART` package v2.9.9 (Sparapani et al. 2021) in R v4.5.1 (R Core Team
2025). We used Claude Opus 5 and Claude Opus 5.5 (Anthropic PBC, 2026) with RStudio
v2025.9.1.401 (Posit Team, 2025) to co-develop the BART pipeline and subsequent counterfactual
analyses. All code was reviewed and edited by MMAB.

A separate model was fitted for every watershed × biotic feature. The training data were the
watershed's low-footprint pixels at which the feature was observed, and the prediction targets
were its high-footprint pixels. Projected easting and northing (EPSG:5072) were included as
predictors so that the trees could absorb spatial structure not explained by the abiotic
covariates; because high-footprint pixels lie within the same watershed as the training pixels,
this is spatial interpolation rather than extrapolation. Predictors that were entirely missing
or invariant among the training pixels were dropped. BART cannot make a prediction for a pixel
with any missing predictor value, so where a predictor was only partially missing (e.g.
phenology metrics over open water) we replaced the missing entries with the median of the
training pixels. Where missingness varied among the training pixels we also added a binary
missingness indicator as a predictor, allowing the trees to split on "was missing" when
missingness was itself informative. Where values were missing only at high-footprint pixels the
indicator would carry no training signal, so the median fill was used alone; these pixels were
predominantly open water and are excluded downstream by the prediction mask (XYZ). Random-number
seeds were derived from the watershed and feature indices, so every fit is reproducible.

Biotic predictors that precede y in the hierarchy entered training at their observed values and
prediction at their backfilled posterior means. Each model therefore conditions on the posterior
mean of the features upstream of it: uncertainty in an upstream feature is not propagated into
the draws of features downstream, and draws are not coupled across features.

Continuous biotic features were fitted with `gbart()` (`type = "wbart"`) using m = 50 trees and
700 posterior draws retained after 300 burn-in iterations, without thinning. We set k = 3
(package default 2), shrinking each tree's contribution more strongly than the default: a
conservative choice for models applied to pixels outside their training set. The tree-structure
prior used the package defaults α = 0.95 and β = 2, so a node splits with probability 0.95 at the
root, 0.24 at depth 1 and 0.11 at depth 2. The prior on σ had 3 degrees of freedom and was
calibrated so that the standard deviation of the log-transformed training response lay at its
90th percentile. Anchoring the prior to this marginal standard deviation, rather than the
package default (the residual standard deviation of a linear regression), places more prior
mass on larger σ and so further discourages overfitting. We enabled the Dirichlet prior on
splitting-variable probabilities (`sparse = TRUE`; Linero 2018) at its package defaults: instead
of choosing split variables uniformly, the sampler learns a probability vector over predictors
and concentrates splits on the informative ones, which suits a predictor set of up to XYZ
variables, most of which are expected to be uninformative for any one biotic feature. Each
continuous predictor was discretised to 100 evenly spaced candidate cutpoints (package default).
Posterior draws were back-transformed as exp(·) − 1 draw by draw, so that posterior means on the
original scale are means of back-transformed draws rather than back-transformed means. Because
the Gaussian model is unbounded on the log scale, a back-transformed draw can fall slightly
below zero; such values were set to zero before re-predicting bird density. At each
high-footprint pixel we retained 100 draws of f(xᵢ), selected at random from the 700; the
adequacy of this subsample is assessed in *Backfilling: prediction*.

Six categorical land-cover layers (XYZ) were fitted with `mbart2()` using probit ensembles
(`type = "pbart"`). This function fits one binary ensemble for each land-cover class present in
the watershed's training pixels, each answering "is this pixel class h, or not?". Each binary
ensemble is fitted by data augmentation (Albert and Chib 1993): a latent normal variable is drawn
for each pixel, truncated to be positive if the pixel belongs to class h and negative otherwise,
and the trees are then updated as in the continuous case with the residual variance fixed at
one. Within each posterior draw, `mbart2()` converts the K ensemble outputs at a pixel into class
scores by normalising exp(f_h) across classes; we averaged these scores over draws and took the
highest-scoring class as the backfilled land cover. Because the ensembles are fitted on a probit
scale but combined by exponential normalisation, the scores sum to one but are not calibrated
class probabilities, so we use them only to rank classes and as a relative measure of
classification uncertainty. Multinomial models used 40 trees per class and k = 3 (giving the
sum of trees a prior standard deviation of 1 on the latent scale), with the same sparsity prior,
and retained 500 draws after 150 burn-in iterations, keeping every tenth iteration (the package
default for binary outcomes). We thinned these chains and not the continuous ones because
latent-variable Gibbs samplers mix more slowly, so consecutive draws are more strongly
autocorrelated. Land-cover layers were backfilled after all continuous features. Each was
predicted from the abiotic predictors, the coordinates and all continuous biotic features, and
no land-cover layer was used as a predictor of any other feature.

For each continuous model we summarised in-sample performance as the root mean square error
(RMSE) and mean absolute error (MAE) of the posterior mean against observed values on the
original scale; the median and 95% credible interval of R², computed as 1 − SSE/SST within each
posterior draw on the original scale; the posterior mean and standard deviation of σ (on the
log(y + 1) scale on which it is estimated); and 95% posterior predictive coverage, the proportion
of observed values inside the 2.5–97.5% interval of posterior predictive draws, generated by
adding N(0, σ²) noise to each draw of f on the log scale before back-transforming. For
out-of-sample performance, a random 10% of each watershed's training pixels was withheld before
fitting (the backfilling models are those fitted to the remaining 90%), and we computed RMSE,
MAE and R² of the posterior mean on the withheld pixels. Because withheld pixels are low-footprint
pixels interspersed with the training pixels, this measures predictive skill within the
low-footprint domain, not transfer to high-footprint pixels, whose abiotic environments differ
(see *Low human footprint: not a proxy for pre-industrialization*); transfer is assessed against
independent data (*Backfill validation: third party data*) and by extrapolation diagnostics
(XYZ). For categorical models we report the accuracy of the highest-scoring class and the mean
Shannon entropy (in nats) of the class scores, in-sample and on the withheld pixels, and
confusion matrices on the withheld pixels (Tables XYZ).

Some watershed × feature combinations did not support a BART fit. For continuous features these
were features that were constant across the training pixels (backfilled with that constant),
that took two or fewer distinct values among the pixels remaining after the holdout split, or
whose watershed held a single high-footprint pixel (both backfilled with the training mean). For
categorical features they were layers with a single class among the training pixels, or among
those remaining after the holdout split, and watersheds with a single high-footprint pixel (all
backfilled with the modal class). These combinations carry no posterior uncertainty and are
excluded from the performance tables.
