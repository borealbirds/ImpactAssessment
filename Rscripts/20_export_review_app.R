# ---
# title: Data bundle for the impactexploreR review app
# author: Mannfred Boehm
# ---

# for impactexploreR, the co-authors' Shiny app (a private sibling repo; APP_DIR redirects,
# default ../impactexploreR): 14C's three tables for every year, 18's change maps as ~1600 px
# WebP, and the Canadian BCRs as SVG paths for the region picker, for the 113 reported species.
# every map gets the same crop box (the union of the maps' non-white extents), so the app can
# swap a map's year in place. runs locally after 14C and 18 for every year

library(magick)
library(sf)

app_dir <- Sys.getenv("APP_DIR", file.path("..", "impactexploreR"))
se_dir  <- "data/derived_data/sector_effects"
map_dir <- "output_figures/change_maps"
years   <- seq(1995, 2020, by = 5)
stopifnot(dir.exists(app_dir))
dir.create(file.path(app_dir, "data"), showWarnings = FALSE)
dir.create(file.path(app_dir, "www", "maps"), recursive = TRUE, showWarnings = FALSE)

# -----------------------
# tables
spp <- trimws(readLines(file.path("Rscripts", "species_113.txt")))
spp <- spp[nzchar(spp) & !startsWith(spp, "#")]
nm  <- read.csv("data/raw_data/species_names.csv")
stopifnot(all(spp %in% nm$species))

# 14D's habitat group of each species (Anna Drake's classes), in 19C's order
groups <- c("conifer" = "Conifer", "mixedwood" = "Mixedwood", "deciduous" = "Deciduous",
            "forest no type" = "Forest, no type", "non-forest" = "Non-forest")
master <- read.csv("data/raw_data/anna_data/master_species_list.csv")
master <- master[master$include, ]
if (!setequal(master$species, spp)) stop("14D's included species differ from species_113.txt")
if (length(odd <- setdiff(master$group, names(groups))))
  stop("14D group(s) not in this script: ", paste(odd, collapse = ", "))
habitat <- setNames(groups[master$group], master$species)[spp]

# 2020's CSVs carry no year suffix, the other years' carry _{year} (14B/14C)
csv_path  <- function(stem, y) file.path(se_dir, paste0(stem, if (y != 2020) paste0("_", y), ".csv"))
read_year <- function(stem) do.call(rbind, lapply(years, function(y) {
  d    <- read.csv(csv_path(stem, y))
  miss <- setdiff(spp, d$species)
  if (length(miss)) stop(basename(csv_path(stem, y)), " lacks ", paste(miss, collapse = ", "))
  cbind(year = y, d[d$species %in% spp, ])
}))
summ <- read_year("results_summary")
secs <- read_year("results_sectors")
bcr  <- read_year("results_bcr")
stopifnot(!anyDuplicated(summ[c("year", "species")]), !anyDuplicated(secs[c("year", "species")]),
          !anyDuplicated(bcr[c("year", "species", "bcr")]))

# species in Table S1's order: 2020 % change, most negative first
s20   <- summ[summ$year == 2020, ]
order <- s20$species[order(s20$change_pct)]

# -----------------------
# BCR polygons -> SVG paths
bcr_names <- c(`3` = "Arctic Plains and Mountains", `4` = "Northwestern Interior Forest",
               `5` = "Northern Pacific Rainforest", `6` = "Boreal Taiga Plains",
               `7` = "Taiga Shield and Hudson Plains", `8` = "Boreal Softwood Shield",
               `9` = "Great Basin", `10` = "Northern Rockies", `11` = "Prairie Potholes",
               `12` = "Boreal Hardwood Transition", `13` = "Lower Great Lakes/St. Lawrence Plain",
               `14` = "Atlantic Northern Forest")
v <- st_read("data/raw_data/Regions/BAM_BCR_NationalModel_Unbuffered.shp", quiet = TRUE)
v <- v[v$country == "can", ]
v$bcr <- paste0("can", v$subUnit)
stopifnot(st_crs(v)$epsg == 5072, all(unique(bcr$bcr) %in% v$bcr))   # 18's map projection
num   <- ifelse(v$subUnit >= 40, v$subUnit %/% 10, v$subUnit)
info  <- data.frame(bcr = v$bcr, label = paste0("BCR ", num, ": ", bcr_names[as.character(num)]))
info  <- info[order(v$subUnit), ]

# one polygon per row; drop islands under 300 km2 (invisible at the picker's ~300 px), keep each
# BCR's largest piece whatever its size
p <- st_make_valid(v[c("bcr")])
st_agr(p) <- "constant"
p <- st_cast(st_cast(st_collection_extract(p, "POLYGON"), "MULTIPOLYGON"), "POLYGON")
p$area <- as.numeric(st_area(p)) / 1e6
p <- p[p$area >= 300 | p$area == ave(p$area, p$bcr, FUN = max), ]
p <- st_simplify(p, preserveTopology = TRUE, dTolerance = 8000)
bb <- st_bbox(p)
ring_d <- function(m) {
  xy <- cbind(round((m[, 1] - bb[["xmin"]]) / 1000), round((bb[["ymax"]] - m[, 2]) / 1000))
  xy <- xy[c(TRUE, rowSums(abs(diff(xy))) > 0), , drop = FALSE]   # drop repeats after rounding
  if (nrow(xy) < 3) return("")
  paste0("M", paste(xy[, 1], xy[, 2], sep = ",", collapse = "L"), "Z")
}
paths <- vapply(info$bcr, function(b) paste(vapply(st_geometry(p[p$bcr == b, ]), function(g)
  paste(vapply(unclass(g), ring_d, ""), collapse = ""), ""), collapse = ""), "")
stopifnot(all(nzchar(paths)))
svg <- list(paths = data.frame(bcr = info$bcr, d = unname(paths)),
            viewbox = sprintf("0 0 %d %d", ceiling((bb[["xmax"]] - bb[["xmin"]]) / 1000),
                              ceiling((bb[["ymax"]] - bb[["ymin"]]) / 1000)))
saveRDS(svg, file.path(app_dir, "data", "bcr_svg.rds"))

# -----------------------
# maps
maps <- expand.grid(species = spp, year = years, stringsAsFactors = FALSE)
maps$src <- file.path(map_dir, sprintf("change_grid_%s_%d.png", maps$species, maps$year))
if (any(!file.exists(maps$src))) stop("missing maps: ", paste(maps$src[!file.exists(maps$src)],
                                                              collapse = ", "))
maps$out <- file.path(app_dir, "www", "maps", sprintf("%s_%d.webp", maps$species, maps$year))

width <- 1800                      # before the crop, which leaves ~1600 px
pad   <- 24
read_small <- function(f) image_resize(image_read(f), paste0(width, "x"))
extent <- function(im) {           # non-white columns (x) and rows (y), 1-based
  g <- as.integer(image_data(im, "gray"))[, , 1] < 248
  c(range(which(colSums(g) > 0)), range(which(rowSums(g) > 0)))
}

# pass 1: every map's extent and size, so that one crop box fits them all
ext <- matrix(NA_integer_, nrow(maps), 6)
for (k in seq_len(nrow(maps))) {
  im <- read_small(maps$src[k])
  ext[k, ] <- c(extent(im), image_info(im)$width, image_info(im)$height)
  rm(im)
  if (k %% 50 == 0) { gc(verbose = FALSE); message("extent ", k, "/", nrow(maps)) }
}
stopifnot(length(unique(ext[, 5])) == 1, length(unique(ext[, 6])) == 1)
x0 <- max(1, min(ext[, 1]) - pad); x1 <- min(ext[1, 5], max(ext[, 2]) + pad)
y0 <- max(1, min(ext[, 3]) - pad); y1 <- min(ext[1, 6], max(ext[, 4]) + pad)
crop <- sprintf("%dx%d+%d+%d", x1 - x0 + 1, y1 - y0 + 1, x0 - 1, y0 - 1)

# pass 2: crop and write
for (k in seq_len(nrow(maps))) {
  im <- image_crop(read_small(maps$src[k]), crop)
  image_write(im, maps$out[k], format = "webp", quality = 82)
  rm(im)
  if (k %% 50 == 0) { gc(verbose = FALSE); message("wrote ", k, "/", nrow(maps)) }
}
stale <- setdiff(list.files(file.path(app_dir, "www", "maps"), full.names = TRUE), maps$out)
if (length(stale)) file.remove(stale)

# -----------------------
# bundle
res <- list(summary = summ, sectors = secs, bcr = bcr,
            names   = setNames(nm$common_name[match(spp, nm$species)], spp),
            habitat = habitat, habitat_groups = unname(groups),
            order   = order, bcr_info = info, years = years,
            map     = list(width = x1 - x0 + 1, height = y1 - y0 + 1),
            results_date = as.Date(max(file.mtime(csv_path("results_summary", 2020)),
                                       file.mtime(csv_path("results_bcr", 2020)))),
            built   = Sys.Date())
saveRDS(res, file.path(app_dir, "data", "results.rds"))

size <- sum(file.size(list.files(file.path(app_dir, c("data", "www/maps")), full.names = TRUE)))
message(sprintf(paste0("%d species x %d years; %d BCR rows in 2020; %d maps (%s px); %d BCR paths;",
                       " bundle %.1f MB -> %s"),
                length(spp), length(years), sum(bcr$year == 2020), nrow(maps), crop, nrow(svg$paths),
                size / 2^20, normalizePath(app_dir)))
