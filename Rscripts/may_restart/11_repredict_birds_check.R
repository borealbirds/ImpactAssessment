# check which species x BCRs have finished from "repredict_birds.R"

# ---- small helpers ----
.exists_ok <- function(f) file.exists(f) && tryCatch({ terra::nlyr(terra::rast(f)) >= 1 }, error = function(e) FALSE)

.done_bcr <- function(species, year, bcr_code) {
  outdir <- file.path(ia_dir, "density_predictions", species, as.character(year))
  obs_fn <- file.path(outdir, sprintf("%s_%s_%s_observed.tif",   species, bcr_code, year))
  bf_fn  <- file.path(outdir, sprintf("%s_%s_%s_backfilled.tif", species, bcr_code, year))
  .exists_ok(obs_fn) && .exists_ok(bf_fn)
}

.model_exists <- function(species, bcr_code) {
  file.exists(file.path(root, "output", "06_bootstraps", species, sprintf("%s_%s.Rdata", species, bcr_code)))
}

# ---- build the TODO grid (species × BCR) ----
species_vec <- c("CAWA", "OSFL")
yr <- 2020

bcr_ref <- bcr_subbasins_ref |> distinct(bcr_label, bcr_code)

todo_grid <-
  tidyr::expand_grid(species = species_vec, bcr_ref) |>
  mutate(
    done  = vapply(seq_len(n()), \(i) .done_bcr(species[i], yr, bcr_code[i]), logical(1)),
    model = vapply(seq_len(n()), \(i) .model_exists(species[i], bcr_code[i]), logical(1))
  ) |>
  filter(!done, model)           # only run where rasters missing and model exists

# quick info
message(sprintf("Already done: %d  |  To do: %d  |  Missing models: %d",
                sum(!todo_grid$model) + sum(!todo_grid$done & !todo_grid$model),
                nrow(todo_grid),
                sum(!todo_grid$model)))

# if nothing to do, bail
if (nrow(todo_grid) == 0) {
  message("Nothing to do — all requested species/BCRs are complete.")
} else {
  
  # group TODOs by species so each worker does one species at a time
  todo_groups <- split(todo_grid$bcr_label, todo_grid$species)
  
  # parallel plan (set your desired worker count)
  plan(sequential)
  
  opts <- furrr::furrr_options(
    seed = TRUE, scheduling = 1,
    packages = c("dplyr","purrr","terra","stringr","tibble","gbm")
  )
  
  # run only on needed BCRs by passing a FILTERED ref into your function
  res_todo <-
    future_map_dfr(
      names(todo_groups),
      \(sp) {
        labs <- todo_groups[[sp]]
        if (length(labs) == 0) return(NULL)
        predict_bird_density_per_bcr_year(sp, yr, bcr_subbasins_ref |> filter(bcr_label %in% labs))
      },
      .options = opts,
      .progress = TRUE
    )
  
  # res_todo has the subbasin summaries for the newly computed BCRs
  # (your mean rasters were written as side-effects inside the function)
}
