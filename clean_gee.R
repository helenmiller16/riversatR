# Clean up the data pulled by GEE and save it into a nice little file

library(data.table)
files <- list.files("data/external/conus/reach_data_sentinel/", 
                    full.names = TRUE)
data_list <- lapply(files, fread)
data <- rbindlist(data_list)
rm(data_list)
data <- data[complete.cases(data)]

# Convert from "digital numbers"
# Harmonized S2 on GEE:
# The assets contain 12 UINT16 spectral bands representing SR scaled by 10000 
SCALE = .0001
data[, `:=`(
  coastal = SCALE*B1, 
  blue = SCALE*B2, 
  green = SCALE*B3, 
  red = SCALE*B4, 
  red_edge1 = SCALE*B5, 
  red_edge2 = SCALE*B6, 
  red_edge3 = SCALE*B7, 
  nir = SCALE*B8, 
  red_edge4 = SCALE*B8A
)]
data[, c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A")] <- NULL
saveRDS(data, "data/external/conus/reflectance.rds")
