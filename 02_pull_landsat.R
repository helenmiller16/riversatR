basin_name <- Sys.getenv("BASIN_NAME")
basin_dir <- file.path("data/external", basin_name)
output_dir <- file.path(basin_dir, "landsat_2015_10_01-2016_10-01")

BANDS   <- c("coastal", "blue", "green", "red", "nir08", "swir16", "swir22")
# Columns to keep. 
# reach_id, date, and bands will always be kept
COLUMNS <- c("valid_pixels", "eo:cloud_cover")

# PARAMETERS FOR CLEANED FILE ----- 
# Remove observations with fewer pixels than this. Keep this pretty low, 
# and potentially use a higher number in downstream analysis. 
# valid pixels > 0 are kept in raw files, this is just for the cleaned file. 
MIN_VALID_PIXELS <- 10 
# Cleaned data file name
CLEANED_FN <- file.path(basin_dir, "landsat_2015_10_01-2016_10-01.csv")




# PULL LANDSAT -----
source("pull_landsat_one_reach.R")
library(future.apply)
plan("multisession", workers = 8) 

if (!dir.exists(output_dir)) dir.create(output_dir)
files <- list.files(output_dir)
reaches_done <- sub("\\.csv", "", files)
water_polygons <- terra::vect("data/external/mississippi_small/water_polygons_buffered.gpkg")
water_polygons <- water_polygons[water_polygons$grwl_pixels > 100, ]
water_polygons <- water_polygons[!as.character(water_polygons$reach_id) %in% as.numeric(reaches_done), ]
# Cannot pass Spat objects to futures because they use external pointers. 
water_polygons <- sf::st_as_sf(water_polygons)
n <- nrow(water_polygons)
out <- list(); length(out) <- n


for (i  in 1:n) {
  geom <- water_polygons[i, ]
  reach_id <- as.character(geom$reach_id)
  out[[i]] <- future( {
    Sys.setenv(
      AWS_REQUEST_PAYER = "requester",
      AWS_DEFAULT_REGION = "us-west-2"
    )
    g <- terra::vect(geom)
    result <- pull_landsat_one_reach(g, 
                                     start_date = as.Date("2015-10-01"), 
                                     end_date = as.Date("2016-10-01"), 
                                     bands = BANDS)
    data.table::fwrite(result, file.path(output_dir, paste0(reach_id, ".csv")))
  }, 
  label = reach_id)
}


### WRITE SUMMARY FILE ----- 

library(data.table)
files <- list.files(output_dir)
reach_ids <- as.numeric(sub("\\.csv", "", files))
data_list <- lapply(file.path(output_dir, files), fread)
. <- sapply(1:length(reach_ids), \(i) data_list[[i]][, reach_id := reach_ids[i]])
rm(.)
data <- rbindlist(data_list, fill = TRUE, use.names = TRUE)
rm(data_list)
# allow NA coastal because it is not included in landsat 7
bands_check <- BANDS[BANDS != "coastal"]
data <- data[complete.cases(data[, ..bands_check])]

# use the image with the greatest pixel count when multiple observations in one day
data$date <- as.Date(data$datetime)
data <- data[valid_pixels > MIN_VALID_PIXELS]
data[, max.pixel.count := max(valid_pixels), c("reach_id", "date")]
data_filtered <- data[valid_pixels == max.pixel.count]
vars <- c("reach_id", "date", BANDS, COLUMNS)
# If there are multiple images with the same number of pixels take the average. 
filtered <- data_filtered[
  , ..vars][
    , c(lapply(.SD, mean, na.rm = TRUE)), by = c("reach_id", "date")]

means <- data[
  , ..vars][
    , c(lapply(.SD, weighted.mean, w = valid_pixels, na.rm = TRUE), valid_pixels_total = sum(valid_pixels)), by = c("reach_id", "date")]
fwrite(filtered, CLEANED_FN)
