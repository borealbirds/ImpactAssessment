# ---
# title:   Follow-up numbers for bart_production_settings.R
# author:  Mannfred Boehm
# created: November 4, 2026
# ---

# Run from the project root after _summary.R, with BPS_WORK (the per-fit files) and BPS_OUT (its
# output directory, holding arrays/: every species' *_2020_coalition_256_arrays.rds from Fir).
# In order: per-pair ratios; sigma drift, chain spread and short chains; BART's share of national
# variance; a rough bound on the error one chain per subbasin leaves out; ESS of a 100-sweep window.

suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr)})
options(width=200, pillar.width=200)
d <- "data/derived_data/rds_files/bart_production_settings"
rd <- function(x) read.csv(file.path(d, paste0(x, ".csv")))
sh <- function(x, t) {cat("\n==", t, "\n"); print(as.data.frame(x), digits=3)}
arms <- rd("arms"); ho <- rd("holdout"); ps <- rd("psis"); fp <- rd("footprint"); sz <- rd("size"); st <- rd("stored_draws"); bu <- rd("burn_in")

# backfill mean: does the tree count move it?
w <- arms |> select(pair,covariate,ntree,fill_mean,obs_mean,contrast,sigma_late,n_fill) |>
  pivot_wider(id_cols=c(pair,covariate,obs_mean,n_fill), names_from=ntree, values_from=c(fill_mean,contrast,sigma_late))
w <- w |> mutate(fill_diff_pct_contrast = 100*(fill_mean_200-fill_mean_50)/abs(contrast_50), sigma_ratio = sigma_late_200/sigma_late_50)
sh(select(w,pair,covariate,n_fill,obs_mean,fill_mean_50,fill_mean_200,contrast_50,fill_diff_pct_contrast,sigma_ratio),"backfill mean, 200 vs 50 trees (late window, 4 chains)")

# pair-wise ratios on holdout
g <- function(tree, win, what) ho |> filter(ntree==tree, window==win) |> select(pair, v=all_of(what)) |> deframe2()
deframe2 <- function(x) setNames(x$v, x$pair)
cmp <- function(a, b, what="rmse") { r <- a/b[names(a)]; c(median=median(r), min=min(r), max=max(r), better=sum(r<1), n=length(r)) }
W1 <- "08B window"; WL <- "late window"; W2 <- "2 chains pooled, 08B length each"; W4 <- "4 chains pooled, 08B length each"
tab <- list()
for (m in c("rmse","crps")) {
  tab[[length(tab)+1]] <- rbind(
   `200 vs 50 trees, 1 chain 08B`      = cmp(g(200,W1,m), g(50,W1,m)),
   `200 vs 50 trees, late window`      = cmp(g(200,WL,m), g(50,WL,m)),
   `50: late vs 08B window`            = cmp(g(50,WL,m), g(50,W1,m)),
   `200: late vs 08B window`           = cmp(g(200,WL,m), g(200,W1,m)),
   `50: 2 chains vs 1`                 = cmp(g(50,W2,m), g(50,W1,m)),
   `50: 4 chains vs 1`                 = cmp(g(50,W4,m), g(50,W1,m)),
   `200: 4 chains vs 1`                = cmp(g(200,W4,m), g(200,W1,m)),
   `200 1 chain vs 50 2 chains (equal)`= cmp(g(200,W1,m), g(50,W2,m)),
   `50 late (1 chain 2000) vs 50 2ch`  = cmp(g(50,WL,m), g(50,W2,m)),
   `200 late vs 50 4 chains`           = cmp(g(200,WL,m), g(50,W4,m)),
   `200 4ch vs 50 4ch`                 = cmp(g(200,W4,m), g(50,W4,m)))
  cat("\n==", m, "ratio (first/second; <1 = first better)\n"); print(round(tab[[length(tab)]],3))
}
sh(ho |> group_by(ntree,window) |> summarise(n=n(), cov=median(coverage90), cov_f=median(coverage90_f_only), .groups="drop"), "coverage")

# PSIS against the holdout
sh(ps |> mutate(psis_vs_holdout=rmse_psis_loo/rmse_holdout, insample_vs_holdout=rmse_insample/rmse_holdout) |>
  group_by(ntree) |> summarise(n=n(), psis_med=median(psis_vs_holdout), psis_min=min(psis_vs_holdout), psis_max=max(psis_vs_holdout), ins_med=median(insample_vs_holdout), k=median(share_k_over_0.7), k_max=max(share_k_over_0.7), sec=median(psis_seconds)), "PSIS vs holdout")
# does PSIS rank 200 vs 50 as the holdout does?
pw <- ps |> select(pair,ntree,rmse_psis_loo,rmse_holdout) |> filter(ntree %in% c(50,200)) |> pivot_wider(names_from=ntree, values_from=c(rmse_psis_loo,rmse_holdout)) |>
  mutate(psis_ratio=rmse_psis_loo_200/rmse_psis_loo_50, hold_ratio=rmse_holdout_200/rmse_holdout_50, agree=(psis_ratio<1)==(hold_ratio<1))
sh(select(pw,pair,psis_ratio,hold_ratio,agree),"200/50 RMSE by PSIS and by holdout")

# footprint
f2 <- fp |> mutate(obs_ratio=rmse_fp_observed/rmse_random, cas_ratio=rmse_fp_cascaded/rmse_random)
sh(f2 |> group_by(ntree) |> summarise(n=n(), obs_ratio=median(obs_ratio), cas_ratio=median(cas_ratio), r2_r=median(r2_random), r2_o=median(r2_fp_observed), r2_c=median(r2_fp_cascaded), cov_r=median(coverage90_random), cov_o=median(coverage90_fp_observed), cov_c=median(coverage90_fp_cascaded)), "footprint vs random, median")
fw <- fp |> select(pair,covariate,ntree,rmse_fp_observed,rmse_fp_cascaded,crps_fp_observed,crps_fp_cascaded) |> pivot_wider(names_from=ntree, values_from=c(rmse_fp_observed,rmse_fp_cascaded,crps_fp_observed,crps_fp_cascaded)) |>
  mutate(obs=rmse_fp_observed_200/rmse_fp_observed_50, cas=rmse_fp_cascaded_200/rmse_fp_cascaded_50, crps_o=crps_fp_observed_200/crps_fp_observed_50, crps_c=crps_fp_cascaded_200/crps_fp_cascaded_50)
sh(select(fw,pair,covariate,obs,cas,crps_o,crps_c),"footprint, 200/50 trees")
cat("median ratios:", round(c(median(fw$obs),median(fw$cas),median(fw$crps_o),median(fw$crps_c)),3), "\n")

sh(sz |> group_by(pair) |> summarise(a=mean(rmse_90),b=mean(rmse_100)) |> mutate(pct=100*(b/a-1)) |> summarise(median=median(pct),mean=mean(pct),min=min(pct),max=max(pct)),"90 to 100% RMSE change")


suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr)})
options(width=200, pillar.width=200)
sh <- function(x, t) {cat("\n==", t, "\n"); print(as.data.frame(x), digits=3)}
wd <- Sys.getenv("BPS_WORK"); nskip <- 300L
crps_sample <- function(x, y) { n <- length(x); xs <- sort(x); mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * xs) / n^2 }
score <- function(f, sigma, y) {
  pred <- colMeans(expm1(f)); noisy <- f + rnorm(length(f), 0, rep(sigma, times=ncol(f))); z <- pmax(expm1(noisy),0)
  lo <- apply(z,2,quantile,.05); hi <- apply(z,2,quantile,.95); flo <- apply(expm1(f),2,quantile,.05); fhi <- apply(expm1(f),2,quantile,.95)
  c(rmse=sqrt(mean((pred-y)^2)), crps=mean(sapply(seq_along(y), function(i) crps_sample(z[,i],y[i]))), cov=mean(y>=lo&y<=hi), cov_f=mean(y>=flo&y<=fhi))
}
ids <- sprintf("pair_%02d",1:8); sig_rows <- list(); chain_rows <- list(); win_rows <- list()
set.seed(2)
for (id in ids) {
  meta <- readRDS(file.path(wd, paste0(id,"_meta.rds"))); nh <- meta$n_fill_sample; y <- meta$y_scored; hf <- seq_len(nh); ho <- nh+seq_along(y)
  om <- mean(meta$obs_fill, na.rm=TRUE)
  for (tr in c(50L,200L)) {
    fs <- list.files(wd, sprintf("^%s_t%03d_c[0-9][.]rds$",id,tr), full.names=TRUE); ch <- map(fs, readRDS)
    K <- min(map_int(ch, function(c) nrow(c$f)))
    # sigma over the chain
    n <- min(lengths(map(ch,"sigma"))); m <- rowMeans(sapply(ch, function(c) c$sigma[1:n]))
    win <- list(`301-800`=301:800,`801-1300`=801:1300,`1301-1650`=1301:1650,`1651-2000`=1651:2000)
    v <- sapply(win, function(r) mean(m[r])); within <- mean(sapply(ch, function(c) sd(c$sigma[1651:2000])))
    sig_rows[[length(sig_rows)+1]] <- c(list(pair=id,ntree=tr), as.list(round(100*(v/v[4]-1),2)), list(last_two_sd=(v[3]-v[4])/within))
    # one chain's backfill mean against another's, in the 08B window
    fm <- map_dbl(ch, function(c) mean(expm1(c$f[1:700,hf])))
    fl <- map_dbl(ch, function(c) mean(expm1(c$f[1001:1700,hf])))
    contrast <- mean(fl)-om
    chain_rows[[length(chain_rows)+1]] <- tibble(pair=id, covariate=meta$pair$covariate, ntree=tr, obs=om, fill=mean(fl), contrast=contrast,
       chain_sd_08B=sd(fm), chain_sd_late=sd(fl), pct_of_fill=100*sd(fm)/mean(fm), pct_of_contrast=100*sd(fm)/abs(contrast), mean_08B=mean(fm), drift_pct_contrast=100*(mean(fm)-mean(fl))/abs(contrast))
    # short chains: four chains pooled, a 100-sweep window after each burn-in
    for (b in c(0L,200L,600L,1000L)) {
      rows <- b+1:100
      fp <- do.call(rbind, map(ch, function(c) c$f[rows,ho,drop=FALSE])); sg <- unlist(map(ch, function(c) c$sigma[nskip+rows]))
      win_rows[[length(win_rows)+1]] <- c(list(pair=id,ntree=tr,burn=nskip+b,len=100L), as.list(score(fp,sg,y)))
    }
    rows <- 1:700; fp <- do.call(rbind, map(ch, function(c) c$f[rows,ho,drop=FALSE])); sg <- unlist(map(ch, function(c) c$sigma[nskip+rows]))
    win_rows[[length(win_rows)+1]] <- c(list(pair=id,ntree=tr,burn=300L,len=700L), as.list(score(fp,sg,y)))
    rm(ch); gc(FALSE)
  }
}
sh(bind_rows(sig_rows),"sigma by window: % above sweeps 1651-2000; last column in within-chain SDs")
cr <- bind_rows(chain_rows); sh(cr,"chain means"); write.csv(cr, file.path(Sys.getenv("BPS_OUT"),"chain_means.csv"), row.names=FALSE)
sh(cr |> group_by(ntree) |> summarise(pct_fill=median(pct_of_fill), pct_contrast=median(pct_of_contrast), drift=median(abs(drift_pct_contrast)), drift_signed=median(drift_pct_contrast)),"medians")
wr <- bind_rows(win_rows); write.csv(wr, file.path(Sys.getenv("BPS_OUT"),"short_chains.csv"), row.names=FALSE)
ref <- wr |> filter(burn==300,len==700) |> select(pair,ntree,rmse0=rmse,crps0=crps)
sh(wr |> left_join(ref,by=c("pair","ntree")) |> group_by(ntree,burn,len) |> summarise(rmse_vs_ref=median(rmse/rmse0), rmin=min(rmse/rmse0), rmax=max(rmse/rmse0), crps_vs_ref=median(crps/crps0), cov=median(cov), cov_f=median(cov_f), .groups="drop"), "four chains pooled, by window (ref = burn 300, 700 kept)")

# BART's share of the reported variance
ad <- file.path(Sys.getenv("BPS_OUT"),"arrays"); fl <- list.files(ad, "coalition_256_arrays", full.names=TRUE)
vs <- map_dfr(fl, function(f) { a <- readRDS(f); ch <- a$obs_on_coal_mat-a$bf_on_coal_mat; cf <- a$obs_total_mat-a$obs_on_coal_mat+a$bf_on_coal_mat
  sp <- sub("_.*","",basename(f))
  w <- function(m) mean(apply(m,1,var)); t <- function(m) var(as.vector(m))
  tibble(species=sp, change_mean=mean(ch), change_sd=sd(as.vector(ch)), bart_share_change=w(ch)/t(ch), bart_share_cf=w(cf)/t(cf), bart_sd_change=sqrt(w(ch))) })
write.csv(vs, file.path(Sys.getenv("BPS_OUT"),"bart_variance_share.csv"), row.names=FALSE)
sh(vs |> summarise(n=n(), q10=quantile(bart_share_change,.1), median=median(bart_share_change), q90=quantile(bart_share_change,.9), max=max(bart_share_change), cf_median=median(bart_share_cf), cf_q90=quantile(bart_share_cf,.9)),"BART (between-scenario) share of the variance of national change and no-industry total")
for (r in c(0.611,0.723)) { infl <- 1/r^2; x <- sqrt(vs$bart_share_change*infl + 1 - vs$bart_share_change); cat(sprintf("sd_ratio %.3f: variance x%.2f; total SD grows median %.1f%%, q90 %.1f%%, max %.1f%%\n", r, infl, 100*(median(x)-1), 100*(quantile(x,.9)-1), 100*(max(x)-1))) }

suppressPackageStartupMessages({library(dplyr);library(readr)})
options(width=200, pillar.width=200)
od <- "data/derived_data/rds_files/bart_production_settings"
v <- read_csv(file.path(od,"bart_variance_share.csv"), show_col_types=FALSE)
s <- read_csv("data/derived_data/sector_effects/shapley_subbasin.csv", col_select=c(species,subbasin,shapley_vegetation_mean,shapley_mean), show_col_types=FALSE)
vs <- s |> group_by(species,subbasin) |> summarise(veg=sum(shapley_vegetation_mean,na.rm=TRUE), tot=sum(shapley_mean,na.rm=TRUE), .groups="drop")
nat <- vs |> group_by(species) |> summarise(veg_national=sum(veg), tot_national=sum(tot), root_sum_sq=sqrt(sum(veg^2)), n_eff=sum(abs(veg))^2/sum(veg^2), .groups="drop")
r <- 0.30
x <- inner_join(v, nat, by="species") |> mutate(miss_sd=r*root_sum_sq, sd_growth=100*(sqrt(1+(miss_sd/change_sd)^2)-1), miss_pct_change=100*miss_sd/abs(tot_national), miss_pct_veg=100*miss_sd/abs(veg_national))
write_csv(x, file.path(od,"national_chain_bound.csv"))
cat("species", nrow(x), "\n")
print(as.data.frame(x |> summarise(across(c(n_eff,sd_growth,miss_pct_change,miss_pct_veg), list(q10=~quantile(.x,.1),med=median,q90=~quantile(.x,.9),max=max)))), digits=3)
cat("\nspecies with SD growth over 10%:", sum(x$sd_growth>10), " over 25%:", sum(x$sd_growth>25), "\n")
cat("sign of national change could flip (|change| < 2 x miss_sd):", sum(abs(x$tot_national)<2*x$miss_sd), "; already within 2 reported SD:", sum(abs(x$tot_national)<2*x$change_sd), "\n")
print(as.data.frame(x |> arrange(desc(bart_share_change)) |> select(species,change_mean,change_sd,bart_share_change,tot_national,veg_national,miss_sd,sd_growth,miss_pct_change) |> head(12)), digits=3)

suppressPackageStartupMessages({library(dplyr);library(purrr)})
options(width=200)
wd <- Sys.getenv("BPS_WORK"); rows <- list()
for (id in sprintf("pair_%02d",1:8)) {
  nh <- readRDS(file.path(wd, paste0(id,"_meta.rds")))$n_fill_sample; hf <- seq_len(nh)
  for (tr in c(50L,200L)) {
    ch <- map(list.files(wd, sprintf("^%s_t%03d_c[0-9][.]rds$",id,tr), full.names=TRUE), function(f) readRDS(f)$f[,hf,drop=FALSE])
    K <- min(map_int(ch,nrow)); pv <- apply(do.call(rbind, map(ch, function(m) m[1001:K,,drop=FALSE])),2,var)
    for (w in list(list("burn 500, 100 kept",201:300), list("burn 300, 700 kept",1:700))) {
      mm <- map(ch, function(m) colMeans(m[w[[2]],,drop=FALSE])); pr <- combn(4,2,simplify=FALSE)
      sq <- Reduce(`+`, map(pr, function(p) (mm[[p[1]]]-mm[[p[2]]])^2))/length(pr)
      sdr <- Reduce(`+`, map(ch, function(m) apply(m[w[[2]],,drop=FALSE],2,sd)))/4/sqrt(pv)
      rows[[length(rows)+1]] <- tibble(pair=id,ntree=tr,window=w[[1]],ess=median(pv/(sq/2)),sd_ratio=median(sdr))
    }
    rm(ch); gc(FALSE)
  }
}
x <- bind_rows(rows); print(as.data.frame(x |> group_by(ntree,window) |> summarise(ess=median(ess),q10=quantile(ess,.1),min=min(ess),sd_ratio=median(sd_ratio),.groups="drop")),digits=3)
