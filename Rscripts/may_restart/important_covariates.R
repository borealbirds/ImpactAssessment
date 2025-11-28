library(BAMexploreR)

data("bam_predictor_importance_v5")

top5_by_species <- 
  bam_predictor_importance_v5 %>%
  group_by(spp, bcr) %>%
  slice_max(order_by = mean_rel_inf, n = 5, with_ties = FALSE) 


predictor_summary <- 
  top5_by_species %>%
  group_by(predictor) %>%
  summarise(
    n_species  = n_distinct(spp),   # number of species it's in the top 5 for
    n_bcrs     = n_distinct(bcr),   # number of BCRs it's associated with
    n_total    = n()                # total number of rows (species–BCR occurrences)
  ) %>%
  arrange(desc(n_total))  

write.csv(predictor_summary, "important_covariates.csv")  
