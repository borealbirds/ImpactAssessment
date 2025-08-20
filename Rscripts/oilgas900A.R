library(sf)
library(tidyverse)

gdb_path <- "C:/Users/mannf/Downloads/900A.gdb"
sf::st_layers(gdb_path)


oilgas_sf <- sf::st_read(gdb_path, layer = "OilAndGas")

