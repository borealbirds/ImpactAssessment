library(tidyverse)
library(sf)


# load data (watch out! big ol' data object about to eat your memory)
root <- "G:/Shared drives/BAM_NationalModels5"
load(file.path(root, "data", "02_NM5.0_data_offsets.Rdata"))

# don't need this much data to demonstrate the bug
rm(bird)
gc()
identical(order(visit$id[1:100]), order(offsets$id[1:100])) #TRUE

# order is matched (as per 02, so safe to subset)
# randomly sample with hopes of catching both CAN and USA surveys
set.seed(777)
visit <- 
  visit |> 
  slice_sample(n=30000)

set.seed(777)
offsets <- 
  offsets |>
  slice_sample(n=30000)

identical(order(visit$id), order(offsets$id)) #TRUE


can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) |> 
  dplyr::select(geometry)


visit.country <- visit |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_intersection(can) |> 
  mutate(country = "can") |> 
  full_join(visit) |> 
  mutate(country = ifelse(is.na(country), "usa", country)) |> 
  st_drop_geometry()

identical(order(visit$id), order(visit.country$id)) #FALSE 

# doesn't match :-(
head(visit$id)
head(visit.country$id)


######
# Anna's example
######

#Test V5 script version

bcr.i="can5"
i=1
spp.i="AMRO"

root2<-"G:/Shared drives/BAM_NationalModels5/data"

load(file.path(root2, "04_NM5.0_data_stratify.Rdata"))

 

covsnew.i <- covlist |>

  dplyr::filter(bcr==bcr.i) |>

  pivot_longer(-bcr, names_to="var", values_to="use") |>

  dplyr::filter(use==TRUE)

 

#6. Get visits to include----

visit.i <- bcrlist[bcrlist[,bcr.i]>0, c("id", bcr.i)] |>

  dplyr::filter(id %in% bootlist[[i + 2]])

 

#7. Get response data (bird data)----

bird.i <- bird[as.character(visit.i$id), spp.i]

 

#8. Get covariates and remove the nonsignificant ones----

cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% c("id",covsnew.i$var)]|> rename(idcov=id)  ### ***checking change:  retain covariate ID ----

 

#9. Add year----

year.i <- visit[visit$id %in% visit.i$id, c("id","year")]|> rename(idyear=id)    ### ****checking change:  retain covariate ID ----

 

#10. Put together data object –

dat.i <- cbind(bird.i, year.i, cov.i) |>

  rename(count = bird.i)|>mutate(birdid=rownames(as.data.frame(bird.i)))  ### ***checking change:  get the bird ID from the matrix rownames----

 

#11. Get offsets----

off.i <- offsets[offsets$id %in% visit.i$id, c("id",spp.i)]|> rename(idoffset=id)
all(off.i$idoffset == dat.i$idcov) 

chk <- cbind(off.i,dat.i)    ### ********checking change:  join the offsets

chk<-chk[c(1:6,61)]  # reduce to the IDs

View(chk)

# but then see: 
table(chk$idoffset == chk$idcov)   # FALSE 24930, TRUE 24
chk[20:30, c("idoffset", "idcov")]