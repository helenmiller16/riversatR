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
saveRDS(data, "data/external/conus/reflectance_raw.rds")

# now cleanup 
data <- readRDS("data/external/conus/reflectance_raw.rds")
# use the image with the greatest pixel count when multiple observations in one day
data$date <- as.Date(data$date)
data <- data[pixel_count > 10]
data[, max.pixel.count := max(pixel_count), c("reach_id", "date")]
data_filtered <- data[pixel_count == max.pixel.count]
filtered <- data_filtered[
  , !c('azimuth', 'zenith', 'tile', 'cloud_cover', 'grwl_pixel_count', 'grwl_code')][
    , c(lapply(.SD, mean)), by = c("reach_id", "date")]
means <- data[
  , !c('azimuth', 'zenith', 'tile', 'cloud_cover', 'grwl_pixel_count', 'grwl_code')][
  , c(lapply(.SD, weighted.mean, w = pixel_count), pixel_count_total = sum(pixel_count)), by = c("reach_id", "date")]
fwrite(filtered, "data/external/conus/reflectance.tsv")
