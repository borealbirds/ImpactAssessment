# ---
# title: Impact Assessment: test if presence of mines is associated with SCANFI "rock" and/or reduced biomass
# author: Mannfred Boehm
# created: June 11, 2025
# ---


root <- "G:/Shared drives/BAM_NationalModels5"


stack_bcr10_2020 <- terra::rast(file.path(root, "gis", "stacks", "can10_2020.tif"))
stack_bcr10_1990 <- terra::rast(file.path(root, "gis", "stacks", "can10_1990.tif"))
mines_vector_mask <- terra::vect(file.path(root, "gis", "other_landscape_covariates", "mincan_dataset_masked.gpkg"))

# extract values from SCANFI a mine locations
mines_vector_mask$scanfi2020 <- terra::extract(stack_bcr10_2020$SCANFI_1km, mines_vector_mask)[,2]
mines_vector_mask$scanfi1990 <- terra::extract(stack_bcr10_1990$SCANFI_1km, mines_vector_mask)[,2]

mines_vector_mask$canhf2020 <- terra::extract(stack_bcr10_2020$CanHF_1km, mines_vector_mask)[,2]
mines_vector_mask$canhf20205x5 <- terra::extract(stack_bcr10_2020$CanHF_5x5, mines_vector_mask)[,2]


# mines_vector_mask$scanfi_cats <- c("bryoid", "herb", "rock", "shrub", "broadleaf", "conifer", "mixed", "water")
table(mines_vector_mask$scanfi2020, useNA = "ifany")
table(mines_vector_mask$scanfi1990, useNA = "ifany")

mines_vector_mask$canhf2020 |> hist(breaks=100)
mean(mines_vector_mask$canhf2020, na.rm=TRUE) #6.82
sd(mines_vector_mask$canhf2020, na.rm=TRUE) #6.66

mines_vector_mask$canhf20205x5 |> hist(breaks=100)
mean(mines_vector_mask$canhf20205x5, na.rm=TRUE) #5.19
sd(mines_vector_mask$canhf20205x5, na.rm=TRUE) #3.72


# find which mines have low HF
mines_df <- as_tibble(mines_vector_mask)

lowhf_mines <- 
  mines_df |>
  filter(canhf20205x5 < 5.2) |>
  select(namemine, canhf20205x5, open1, close1, open2, close2, open3, close3)


mine_hf_3x3 <- 
  terra::extract(stack_bcr10_2020$CanHF_1km, 
                              mines_vector_mask, 
                              buffer = 1500,  # 1.5 km radius should capture 3x3 for 1 km resolution
                              fun = NULL,     # return raw pixel values
                              cells = FALSE)
mine_hf_summary <- mine_hf_3x3 |>
  group_by(ID) |>
  summarize(hf_min = min(CanHF_1km, na.rm = TRUE),
            hf_max = max(CanHF_1km, na.rm = TRUE),
            hf_mean = mean(CanHF_1km, na.rm = TRUE))

