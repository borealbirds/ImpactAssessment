library(tidyverse)
library(sf)
library(rgee)
ee_Initialize() # or however you get rgee going on your local

# load data 
root <- "G:/Shared drives/BAM_NationalModels5"

# 03's loc.n: one row per location x year, the surveys we'll get values for from GEE
loc.n <- data.table::fread(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD.csv"),
                           select=c("project", "location", "lat", "lon", "year", "id"), data.table=FALSE)

# not in year order
is.unsorted(loc.n$year)

# don't need this much data to demonstrate the bug
set.seed(777)
loc.n <-
  loc.n |>
  slice_sample(n=100) |>
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |>
  st_transform(crs=5072)


# 03:651-652, how 03 starts loc.gee before any layer is saved: a plain copy of loc.n
loc.gee <- data.frame(loc.n) |>
  dplyr::select(-geometry)

# 03:52 and 03:648, the layer lookup and filter for year-matched covariates
meth <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", TemporalResolution=="match")

# one pass of 03's layer loop (03:657), for summer temperature
i <- which(meth.gee$Label=="ERATavesm_1km")

# this portion is copy-pasted from 03_ExtractCovariates.R (line 660-669)
years.gee <- seq(meth.gee$GEEYearMin[i], meth.gee$GEEYearMax[i])
dt = data.table::data.table(year=years.gee, val=years.gee)
data.table::setattr(dt, "sorted", "year")
data.table::setkey(dt, year)
loc.n.i <- loc.n
loc.n.i$yearrd <- dt[J(loc.n$year), roll = "nearest"]$val

# lines 673-750: extract one year at a time, in batches of 100
# loop takes ~20 seconds
# some elements are NULL because our 100 sample survey doesn't cover every possible year
loc.j <- list()
for(j in 1:length(years.gee)){

  loc.n.yr <- dplyr::filter(loc.n.i, yearrd==years.gee[j]) |>
    mutate(loop = ceiling(row_number()/100))

  if(nrow(loc.n.yr) > 0){

    loc.k <- list()
    for(k in 1:max(loc.n.yr$loop)){

      loc.n.loop <- dplyr::filter(loc.n.yr, loop==k)
      point <- sf_as_ee(loc.n.loop)

      start.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-", meth.gee$GEEMonthMin[i], "-01")

      if(meth.gee$GEEMonthMax[i] > meth.gee$GEEMonthMin[i]){
        end.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-", meth.gee$GEEMonthMax[i], "-28")
      }
      if(meth.gee$GEEMonthMax[i] < meth.gee$GEEMonthMin[i]){
        end.k <- paste0(years.gee[j], "-", meth.gee$GEEMonthMax[i], "-28")
      }

      img.i <- ee$ImageCollection(meth.gee$Link[i])$filter(ee$Filter$date(start.k, end.k))$select(meth.gee$GEEBand[i])$mean()
      loc.k[[k]] <- ee_extract(x=img.i, y=point, scale=meth$GEEScale[i])
    }

    loc.j[[j]] <- data.table::rbindlist(loc.k, fill=TRUE)
    colnames(loc.j[[j]]) <- c(colnames(loc.j[[j]])[1:ncol(loc.j[[j]])-1], meth.gee$Label[i])
  }
}

# lines 769-770, stack the years, then bind by position
loc.i <- data.table::rbindlist(loc.j, fill=TRUE)
loc.gee <- cbind(loc.gee, loc.i |>
                   dplyr::select(meth.gee$Label[i]))

# row order does not match
identical(loc.gee$id, loc.i$id) #FALSE

# doesn't match :-(
head(loc.gee$id)
head(loc.i$id)

# each survey's own value, joined by id
# `own` is each survey's correctly extracted average temperature 
loc.gee <- left_join(loc.gee, loc.i |> dplyr::select(id, own=ERATavesm_1km))

# if the cbind kept each value on its own survey, this would be 1
cor(loc.gee$ERATavesm_1km, loc.gee$own) # 0.016, almost no correlation
