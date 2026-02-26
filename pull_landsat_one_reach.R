# TODO: set up gdal config file
# https://gdal.org/en/stable/user/configoptions.html#gdal-configuration-file

# Pull Landsat data for one reach given an NHD id, GRWL id, SWOT id, or lat/lon
# Currently only coordinates supported. 

pull_landsat_one_reach <- function(
    # nhd_id = NULL, 
    # grwl_id = NULL, 
    # SWOT_id = NULL, 
    coords = NULL, # c(lon, lat)
    reach_polygons_file = NULL,
    start_date, 
    end_date, 
    tolerance = 20) {
  
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
  .landsat_url <- function(x) {
    sub("s3://", "/vsis3/", x)
  }
  
  # Get matching reach polygon
  reach_polygons <- sf::read_sf(reach_polygons_file)
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
    rstac::items_filter(properties$`eo:cloud_cover` < 10)
  
  # set up data frame to store results
  assets_of_interest <- c("blue", "green","red","nir08","swir16", "swir22","qa_pixel")
  data_list <- lapply(1:items_length(items), \(i) {
    item <- rstac::items_select(items, i)
    urls <- rstac::assets_url(item, assets_of_interest)
    r <- terra::rast(landsat_url(urls))
    names(r) <- assets_of_interest
    r <- crop(r, sf::st_transform(polygon, crs(r)))
    # remove cloudy values
    qa_mask <- rsi::landsat_mask_function(r$qa_pixel, include = "water")
    # count the number of valid pixels
    valid_pixels <- terra::zonal(qa_mask, vect(polygon), fun = sum, na.rm = TRUE)
    names(valid_pixels) <- "valid_pixels"
    r <- terra::mask(r, qa_mask, maskvalues = FALSE)
    band_medians <- terra::zonal(r, terra::vect(sf::st_transform(polygon, terra::crs(r))), fun = median, na.rm = TRUE)
    band_medians["qa_pixel"] <- NULL
    band_medians <- band_medians * 0.0000275 - 0.2
    cbind(valid_pixels, band_medians)
  })
  Sys.setenv(
    AWS_REQUEST_PAYER = AWS_REQUEST_PAYER,
    AWS_DEFAULT_REGION = AWS_DEFAULT_REGION
  )
  data_dt <- data.table::rbindlist(data_list)
  items_dt <- data.table::as.data.table(rstac::items_as_tibble(items))
  data <- data.table::cbind(items_dt, data_dt)
  data$`processing:software` <- sapply(data$`processing:software`, \(x) paste(x, collapse = "_"))
  data$`proj:centroid` <- sapply(data$`proj:centroid`, \(x) paste(x, collapse = " "))
  data
}
