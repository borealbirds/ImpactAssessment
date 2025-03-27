# ---
# title: National Models 5.0 - predict bird abundance in backfilled landscape
# author: Mannfred Boehm
# created: March 25, 2025
# ---

#1. preamble----

library(tidyverse)
library(terra)
library(gbm)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"



#2. import raster stacks----

stack_bcr14_2020_bf <- terra::rast(file.path(getwd(), "data/derived_data/backfilled_rasters/can14_2020_backfilled.tif"))
stack_bcr14_2020 <- terra::rast(file.path(root, "gis", "stacks", "can14_2020.tif"))


#3. import bird density model----

my_colours <- colorRampPalette(c( "#F0E442", "#009E73", "#0072B2"))(10)


# loads gbm model as `b.i`
load(file=file.path(root, "output", "bootstraps", "CAWA", "CAWA_can14_1.R"))

prebackfill_pop  <- terra::predict(model=b.i, object=stack_bcr14_2020, type="response")
postbackfill_pop <- terra::predict(model=b.i, object=stack_bcr14_2020_bf, type="response")

# mask predictions and count birds
prebackfill_pop2 <- terra::mask(prebackfill_pop, stack_bcr14_2020$year)


postbackfill_pop2 <- terra::mask(postbackfill_pop, stack_bcr14_2020$year)


# check distribution of values and plot
quantile(values(prebackfill_pop2), na.rm=TRUE)
threshold <- as.numeric(quantile(values(prebackfill_pop2), probs = 0.8, na.rm = TRUE))
terra::values(prebackfill_pop2) <- pmin(values(prebackfill_pop2), threshold)
sum(terra::values(prebackfill_pop2)*100, na.rm=TRUE) #653359.6

plot(prebackfill_pop2, col=my_colours)



quantile(values(postbackfill_pop2), na.rm=TRUE)
# threshold <- as.numeric(quantile(values(postbackfill_pop2), probs = 0.8, na.rm = TRUE))
terra::values(postbackfill_pop2) <- pmin(values(postbackfill_pop2), threshold)
sum(terra::values(postbackfill_pop2)*100, na.rm=TRUE) #665597.7 (+12,238)

plot(postbackfill_pop2, col=my_colours)



# plotting some covariates
breaks <- seq(0, 160, by = 40)
plot(stack_bcr14_2020$SCANFIbiomass_1km, breaks=breaks)
plot(stack_bcr14_2020_bf$SCANFIbiomass_1km, breaks=breaks)
