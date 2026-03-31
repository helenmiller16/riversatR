terraOptions(memfrac=.8)


## STEP 2: Pull data at reaches associated with my sensor locations. 

# TODO: set up gdal config file
# https://gdal.org/en/stable/user/configoptions.html#gdal-configuration-file

# Pull Landsat data for one reach given an NHD id, GRWL id, SWOT id, or lat/lon
# Currently only coordinates supported. 

.landsat_url <- function(x) {
  sub("s3://", "/vsis3/", x)
}

# include = water, land, both
.landsat_qa_mask <- function(qa,
                             include = "water") {
  # Make masks
  dilated_cloud <- 2 # bit 1
  cirrus <- 4 # bit 2
  cloud <- 8 # bit 3
  cloud_shadow <- 16 # bit 4
  snow <- 32 # bit 5
  clear <- 64 # bit 6
  water <- 128 # bit 7
  
  q <- terra::values(qa)
  
  m <- bitwAnd(q, cloud) == 0 & 
    bitwAnd(q, cloud_shadow) == 0 & 
    bitwAnd(q, dilated_cloud) == 0 & 
    bitwAnd(q, cirrus) == 0 & 
    bitwAnd(q, snow) == 0
  if (include == "water") {
    m <- bitwAnd(q, water) != 0
  } else if (include == "land") {
    m <- bitwAnd(q, water) == 0
  }
  mask <- terra::rast(qa)
  terra::values(mask) <- m
  mask
}

pull_landsat_one_reach <- function(
    # nhd_id = NULL, 
  # grwl_id = NULL, 
  # SWOT_id = NULL, 
  coords = NULL, # c(lon, lat)
  reach_polygons_file = NULL,
  start_date, 
  end_date, 
  tolerance = 20, 
  # max_cloud_cover = 100, # hardcoded because this didn't work in the function for some reason
  bands = c("coastal", "blue", "green","red","nir08","swir16", "swir22")) {
  
  # reach_specified <- sum(c(!is.null(nhd_id), !is.null(grwl_id), !is.null(SWOT_id), !is.null(coords)))
  # if (reach_specified == 0) {
  #   stop("specify a reach.")
  # } else if (reach_specified > 1) {
  #   stop("specify a reach only one way.")
  # }
  # 
  # # If NHDplusHR reach id: 
  # if (!is.null(nhd_id)) {
  #   library(nhdplusTools)
  #   huc <- get_huc()
  # } 
  units(tolerance) <- "m"
  # coords <- c(-119.359875, 46.582473)
  # reach_polygons_file <- "data/external/columbia/water_polygons.gpkg"
  # start_date <- as.Date("2011-01-01")
  # end_date <- as.Date("2011-04-15")
  if (!is.null(coords)) {
    location <- sf::st_sfc(st_point(coords), crs = "EPSG:4326")
  } else {
    stop("specify coords")
  }
  
  AWS_REQUEST_PAYER <- Sys.getenv("AWS_REQUEST_PAYER")
  AWS_DEFAULT_REGION <- Sys.getenv("AWS_DEFAULT_REGION")
  Sys.setenv(
    AWS_REQUEST_PAYER = "requester",
    AWS_DEFAULT_REGION = "us-west-2"
  )
  
  # Get matching reach polygon
  reach_polygons$id <- 1:nrow(reach_polygons)
  location <- sf::st_transform(location, sf::st_crs(reach_polygons))
  # get closest reach
  polygon <- reach_polygons[sf::st_nearest_feature(location, reach_polygons), ]
  if (sf::st_distance(location, polygon)[[1]] > tolerance) {
    stop(paste0("location is not within ", tolerance, "m of a valid reach."))
  }
  
  date_range <- paste0(format(start_date, "%Y-%m-%dT00:00:00Z"), 
                       "/", 
                       format(end_date, "%Y-%m-%dT00:00:00Z"))
  # pull landsat for that polygon
  stac_client <- rstac::stac("https://earth-search.aws.element84.com/v1")
  collection <- "landsat-c2-l2"
  items <- rstac::stac_search(stac_client,
                              collections = collection,
                              bbox = sf::st_bbox(sf::st_transform(polygon, "EPSG:4326")),
                              datetime = date_range) |>
    rstac::get_request() |>
    rstac::items_fetch() |> 
    rstac::items_filter(properties$`platform` %in% c("landsat-7", "landsat-8", "landsat-9")) |> 
    rstac::items_filter(properties$`eo:cloud_cover` < 90)
  
  data_list <- pbapply::pblapply(1:rstac::items_length(items), \(i) {
    
    item <- rstac::items_select(items, i)
    b <- intersect(bands, rstac::items_assets(item))
    # set up data frame to store results
    assets_of_interest <- c(b, "qa_pixel")
    urls <- rstac::assets_url(item, assets_of_interest)
    r <- terra::rast(.landsat_url(urls))
    names(r) <- assets_of_interest
    if (terra::crs(polygon) != terra::crs(r)) {
      polygon <- sf::st_transform(polygon, terra::crs(r))
    }
    r <- terra::crop(r, polygon)
    # remove cloudy values
    qa_mask <- .landsat_qa_mask(r$qa_pixel, include = "water")
    # also exclude cells adjacent to masked-out cells
    qa_mask <- terra::focal(qa_mask, matrix(1,3,3), all)
    # count the number of valid pixels
    valid_pixels <- terra::zonal(qa_mask, terra::vect(polygon), fun = sum, na.rm = TRUE)
    names(valid_pixels) <- "valid_pixels"
    if (valid_pixels > 0) {
      r <- terra::mask(r, qa_mask, maskvalues = 0)
      band_medians <- terra::zonal(r, terra::vect(sf::st_transform(polygon, terra::crs(r))), fun = median, na.rm = TRUE)
      band_medians["qa_pixel"] <- NULL
      band_medians <- band_medians * 0.0000275 - 0.2
      cbind(valid_pixels, band_medians)
    } else {
      valid_pixels
    }
    
  })
  Sys.setenv(
    AWS_REQUEST_PAYER = AWS_REQUEST_PAYER,
    AWS_DEFAULT_REGION = AWS_DEFAULT_REGION
  )
  data_dt <- data.table::rbindlist(data_list, fill = TRUE)
  items_dt <- data.table::as.data.table(rstac::items_as_tibble(items))
  data <- cbind(items_dt, data_dt)
  data$`processing:software` <- sapply(data$`processing:software`, \(x) paste(x, collapse = "_"))
  data$`proj:centroid` <- sapply(data$`proj:centroid`, \(x) paste(x, collapse = " "))
  data
}


# Plotting helpers
# terra::plotRGB(r * .0000275 - .2, 3,2,1, stretch = "lin", zlim = c(0, .2))


# Plotting helpers
# 4,3,2 for L8 (with coastal band)
# 3,2,1 for L7. 
#terra::plotRGB(terra::subset(r, b)* .0000275 - .2, 4,3,2, stretch = "lin", zlim = c(0, .2))


# Get data matching metabolism site

library(sf)
library(data.table)
source("pull_landsat_one_reach.R")
reach_polygons_file <- "C:/github/river-sat/data/external/conus/water_polygons.gpkg"
reach_polygons <- read_sf(reach_polygons_file)
metab_data <- read_sf("C:/github/metabolism-mississippi/data/external/appling/sites_with_dates.gpkg")

data_list <- lapply(1:nrow(metab_data), \(row) {
  print(row)
  d <- metab_data[row, ]
  coords <- st_coordinates(d)[1, ]
  start_date <- d$start_date
  end_date <- d$end_date
  bands <- c("coastal", "blue", "green","red","nir08","swir16", "swir22")
  try(pull_landsat_one_reach(
    coords = coords, 
    reach_polygons = reach_polygons,
    start_date = start_date, 
    end_date = end_date, 
    tolerance = 100, 
    bands = bands))
})
saveRDS(data_list, "C:/github/river-sat/inst/data_list.rds")
names(data_list) <- metab_data$site_name
data_list2 <- lapply(names(data_list), \(n) {
  d <- data_list[[n]]
  if ("try-error" %in% class(d)) return(NULL)
  d$site_name <- n
  d
})
data <- rbindlist(data_list2)
fwrite(data, "C:/github/metabolism-mississippi/data/external/reflectance/landsat_appling_sites_better.tsv")
# m <- metab_data[grepl("baton", metab_data$long_name, ignore.case = TRUE), ]
# coords <-  st_coordinates(m)[1, ]
# start_date <- m$start_date
# end_date <- m$end_date
# 
# br_data <- pull_landsat_one_reach(
#   coords, 
#   reach_polygons_file = "C:/github/river-sat/data/external/conus/water_polygons.gpkg",
#   start_date, 
#   end_date)
# saveRDS(br_data, "C:/github/river-sat/inst/br_data.rds")
