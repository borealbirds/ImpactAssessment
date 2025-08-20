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

# import BCR14 boundary 
bcr14_boundary <- 
  terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  terra::project(x=_, y=stack_bcr14_2020) |> 
  terra::crop(x=_, y=stack_bcr14_2020)


#3. import bird density model----

# loads gbm model as `b.i`
load(file=file.path(root, "output", "bootstraps", "CAWA", "CAWA_can14_1.R"))

prebackfill_pop  <- terra::predict(model=b.i, object=stack_bcr14_2020, type="response")
postbackfill_pop <- terra::predict(model=b.i, object=stack_bcr14_2020_bf, type="response")

# mask predictions to exclude ocean
prebackfill_pop2 <- terra::mask(prebackfill_pop, stack_bcr14_2020$year)
postbackfill_pop2 <- terra::mask(postbackfill_pop, stack_bcr14_2020$year)



#4. postbackfill analysis---- 


# use median absolute deviation (MAD) to
# remove outliers and estimate population size
med_post <- median(values(postbackfill_pop2), na.rm = TRUE)
mad_post <- mad(values(postbackfill_pop2), na.rm = TRUE)
threshold_post <- med_post + (2.5 * mad_post)

vals_post_clipped <- ifelse(values(postbackfill_pop2) > threshold_post, 0, values(postbackfill_pop2))
terra::values(postbackfill_pop2) <- vals_post_clipped

sum(terra::values(postbackfill_pop2)*100, na.rm=TRUE) #658724 (+19413)
hist(values(postbackfill_pop2),breaks=100,ylim=c(0,20000), yaxp=c(0,20000,5000))





#5. prebackfill analysis----

# use median absolute deviation (MAD) to
# remove outliers and estimate population size
# we'll use the threshold from the post-backfill pop. estimates because
# (for now) the pre-backfill estimates have some weird outliers
vals_pre_clipped <- ifelse(values(prebackfill_pop2) > threshold_post, 0, values(prebackfill_pop2))
terra::values(prebackfill_pop2) <- vals_pre_clipped

sum(terra::values(prebackfill_pop2)*100, na.rm=TRUE) #633406
hist(values(prebackfill_pop2),breaks=100,ylim=c(0,20000), yaxp=c(0,20000,5000))




#6. plotting CAWA distributions----

my_colours <- colorRampPalette(c( "#F0E442", "#009E73", "#0072B2"))(10)
range <- range(values(postbackfill_pop2,na.rm=TRUE))
plot(postbackfill_pop2, col=my_colours, range=range)

plot(prebackfill_pop2, col=my_colours)
lines(bcr14_boundary, col="black", lwd=1)



#7. create a change raster and plot----

my_colours <- colorRampPalette(c( "#0072B2","#EEE9E9","#D55E00"))(5)


change_raster <- postbackfill_pop2 - prebackfill_pop2
plot(change_raster, col=my_colours)
lines(bcr14_boundary, col="black", lwd=1)




#8. miscellaneous plotting for FY24-25 report----

my_colours <- colorRampPalette(c("#F0E442", "#009E73", "#0072B2"))(10)

# canopy height
# set range to the backfilled range (wider than pre-backfilled range)
range <- range(values(stack_bcr14_2020_bf$SCANFIheight_1km,na.rm=TRUE))
#plot(project(stack_bcr14_2020$SCANFIheight_1km,"EPSG:4326"), col=my_colours, range = range)
plot(stack_bcr14_2020$SCANFIheight_1km, col=my_colours, range = range)
plot(stack_bcr14_2020_bf$SCANFIheight_1km, col=my_colours, range = range)
lines(bcr14_boundary, col="black", lwd=1)

# SCANFI
levels(stack_bcr14_2020$SCANFI_1km) <- scanfi_cats
plot(stack_bcr14_2020$SCANFI_1km, col=my_colours)
plot(stack_bcr14_2020_bf$SCANFI_1km, col=my_colours)
lines(bcr14_boundary, col="black", lwd=1)
