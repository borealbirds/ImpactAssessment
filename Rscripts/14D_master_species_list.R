# ---
# title: Master species list with habitat groups
# author: Mannfred Boehm
# created: October 22, 2026
# ---
# One row per species in species_todo.txt: the habitat group its results are reported under,
# and whether it is reported at all (`include`). Groups, in order of precedence:
#   withheld                BAM's release filter keeps none of its BCRs (absent from 14C)
#   wetland, grassland      Birds Canada's lists (UPSA, on both, is grassland)
#   open                    Anna's tundra / above-treeline list
#   non-forest, conifer, mixedwood, deciduous, forest no type   from Anna's classes (below)
# Only the last five are included (113 species on 2026-10-22). The first four lists, and Birds
# Canada's a-priori non-forest list, came in Anna's email of 2026-10-19 and are kept in
# anna_data/species_lists.csv, gitignored like her three class tables.
#
# Anna's classes are per species x BCR x scale (Loc* = local, Land* = landscape). A label counts
# when it comes from one preference window and S > ln(1.5), as in her group analysis. ONE label
# per species, for national and BCR results alike: each BCR's label votes, weighted by the
# species' observed 2020 birds there (14C's results_bcr.csv), in this order of scales:
#   1. local (Anna, 2026-10-22: fewer species had landscape values contribute to model fit)
#   2. landscape
#   3. Birds Canada's non-forest list (forest if the species is not on it)
# Her tables have no can3, can11 or can42. can11 votes with the species' can61 label (Anna: its
# forest species sit along the border with can61); can3 and can42 do not vote. A forest species
# then gets a conifer type (conifer, mixedwood, deciduous) by the same vote, or "forest no type".
#
# Run locally after 14C (anna_data/ is not on the cluster):
#   Rscript Rscripts/14D_master_species_list.R
# Writes data/raw_data/anna_data/master_species_list.csv.
# ---

suppressPackageStartupMessages(library(tidyverse))

ia_dir   <- getwd()   # the local RProject
anna_dir <- file.path(ia_dir, "data", "raw_data", "anna_data")

S_MIN      <- log(1.5)
CAN11_FROM <- "can61"

read_anna <- function(f) read_csv(file.path(anna_dir, f), show_col_types = FALSE,
                                  name_repair = "unique_quiet")

# every species, including the six commented out as "# CAWA (done on the C2d backfill, ...)"
species <- sub("^#?\\s*(\\S*).*$", "\\1", readLines(file.path(ia_dir, "Rscripts", "species_todo.txt")))
species <- unique(species[nzchar(species)])
stopifnot(all(grepl("^[A-Z]{4}$", species)))

lists <- read_anna("species_lists.csv")
on_list <- function(l) lists$species[lists$list == l]

# observed 2020 birds per species x BCR, over the BCRs BAM releases: the vote weights
pop <- read_csv(file.path(ia_dir, "data", "derived_data", "sector_effects", "results_bcr.csv"),
                show_col_types = FALSE) |>
  select(species, bcr, observed)

# usable labels, one row per species x BCR x scale --------------------------------
closure <- read_anna("OpenHabitatIndicators.csv") |>
  filter(type %in% c("Loc_closure", "Land_closure")) |>
  group_by(spp, BCR, type) |> filter(n() == 1L) |> ungroup() |>     # one preference window
  filter(open %in% c("forest", "nonforest"), S > S_MIN) |>
  transmute(species = spp, bcr = BCR, scale = if_else(startsWith(type, "Loc"), "local", "landscape"),
            label = open)

ctype <- read_anna("ForestConifer_Indicators.csv") |>
  filter(type %in% c("Loc_pConifer", "Land_pConifer")) |>
  group_by(spp, BCR, type) |> filter(n() == 1L) |> ungroup() |>
  filter(class %in% c("conifer", "mixedwood", "deciduous"), S > S_MIN) |>
  transmute(species = spp, bcr = BCR, scale = if_else(startsWith(type, "Loc"), "local", "landscape"),
            label = class)

stopifnot(!"can11" %in% c(closure$bcr, ctype$bcr))

# population-weighted vote at one scale: the winning label, its share of the voting birds and
# the number of BCRs that voted
vote <- function(labels, sc) {
  l <- filter(labels, scale == sc) |> select(species, bcr, label)
  l <- bind_rows(l, filter(l, bcr == CAN11_FROM) |> mutate(bcr = "can11"))
  inner_join(pop, l, by = c("species", "bcr")) |>
    group_by(species, label) |>
    summarise(w = sum(observed), n_bcr = n(), .groups = "drop_last") |>
    mutate(share = w / sum(w), n_bcr = sum(n_bcr)) |>
    slice_max(w, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(species, label, share = round(share, 3), n_bcr, source = sc)
}

# the first scale with a vote decides
decide <- function(labels) {
  bind_rows(vote(labels, "local"), vote(labels, "landscape")) |>
    mutate(source = factor(source, c("local", "landscape"))) |>
    group_by(species) |> slice_min(source, n = 1) |> ungroup() |>
    mutate(source = as.character(source))
}
forest <- decide(closure)
ftype  <- decide(ctype)

# assemble ---------------------------------------------------------------------------
reported <- unique(pop$species)

out <- tibble(species = species) |>
  left_join(rename_with(forest, ~ paste0("forest_", .x), -species), by = "species") |>
  left_join(rename_with(ftype,  ~ paste0("type_", .x),   -species), by = "species") |>
  mutate(
    birdscanada   = if_else(species %in% on_list("birdscanada_nonforest"), "nonforest", "forest"),
    forest_source = coalesce(forest_source, "birdscanada"),
    forest        = coalesce(forest_label, birdscanada),
    disagrees_with_birdscanada = forest_source != "birdscanada" & forest != birdscanada,
    type          = if_else(forest == "forest", type_label, NA_character_),
    across(starts_with("type_"), ~ if_else(forest == "forest", .x, NA)),
    group = case_when(
      !species %in% reported            ~ "withheld",
      species %in% on_list("grassland") ~ "grassland",
      species %in% on_list("wetland")   ~ "wetland",
      species %in% on_list("open")      ~ "open",
      forest == "nonforest"             ~ "non-forest",
      is.na(type)                       ~ "forest no type",
      TRUE                              ~ type),
    include = group %in% c("non-forest", "conifer", "mixedwood", "deciduous", "forest no type")) |>
  left_join(summarise(pop, observed_2020 = sum(observed), .by = species), by = "species") |>
  select(species, include, group, forest, forest_source, forest_share, forest_n_bcr,
         type, type_source, type_share, type_n_bcr, birdscanada, disagrees_with_birdscanada,
         observed_2020) |>
  arrange(desc(include), group, species)

stopifnot(!anyDuplicated(out$species), !anyNA(out$group))
missing_list <- setdiff(lists$species, species)
if (length(missing_list) > 0L)
  message("on Anna's lists but not in species_todo.txt (ignored): ", paste(sort(missing_list), collapse = " "))

write_csv(out, file.path(anna_dir, "master_species_list.csv"), na = "")

message("wrote ", file.path(anna_dir, "master_species_list.csv"), ": ", nrow(out), " species, ",
        sum(out$include), " included")
print(count(out, include, group), n = Inf)
message("forest / non-forest decided by: ",
        paste(names(table(out$forest_source[out$include])), table(out$forest_source[out$include]),
              sep = " ", collapse = ", "),
        " | disagree with Birds Canada: ", sum(out$disagrees_with_birdscanada & out$include))
