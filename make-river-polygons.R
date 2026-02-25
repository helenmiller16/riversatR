basin_name <- Sys.getenv("BASIN_NAME")
basin_dir <- file.path("data/external", basin_name)
if (!dir.exists(basin_dir)) dir.create(basin_dir)

crs <- paste0("EPSG:", Sys.getenv("CRS"))
basin_path <- Sys.getenv("WATERSHED_FILE")
# Data prerequisites: 
# Global river widths from Landsat 
# https://zenodo.org/records/1297434
grwl_mask_path <- file.path(Sys.getenv("DATA_DIR"), "GRWL_mask_V01.01")
grwl_tile_path <- file.path(Sys.getenv("DATA_DIR"), "GRWL_tiles")
grwl_basin_vrt <- file.path(basin_dir, paste0(basin_name, "_grwl.vrt"))

# SWORD (choose appropriate continent)
# 10km reaches >30m wide. 
# https://zenodo.org/records/10013982
sword_path <- file.path(Sys.getenv("DATA_DIR"), 
                        "SWORD_v16_gpkg", 
                        Sys.getenv("SWORD_FILE"))

library(sf)
library(terra)
library(matrixStats)

# Get reaches within the river basin
basin <- read_sf(basin_path)
reaches <- read_sf(sword_path)
if (!st_crs(basin) == st_crs(reaches)) {
  basin <- st_transform(basin, st_crs(reaches))
}
reaches <- st_intersection(reaches, basin)
write_sf(reaches, 
         file.path("data/external", basin_name, "reaches.gpkg"))


if (!file.exists(grwl_basin_vrt)) {
    # Use GRWL mask to create river polygons
    # First need to figure out which tiles to use
    grwl_tiles <- read_sf(grwl_tile_path)
    basin_tiles <- st_intersection(grwl_tiles, basin)
    tile_ids <- basin_tiles$TILE_ID
    files <- file.path(grwl_mask_path, paste0(tile_ids, ".tif"))
    # Tiles without any rivers are excluded, I think (e.g. NI10)
    files <- files[file.exists(files)]
    rast_list <- lapply(files, rast)
    # project them all
    # This takes a few minutes but should make everything else 
    # much quicker and more straightforward.
    rast_list <- lapply(rast_list, project, crs, method = "near")
    rast_sprc <- sprc(rast_list)
    # virtual raster dataset for tiles
    grwl_mask <- vrt(rast_sprc, filename = grwl_basin_vrt)
} else {
    grwl_mask <- rast(grwl_basin_vrt)
}


# # buffer each reach based on mean width

reaches <- st_transform(reaches, crs)
geoms <- st_geometry(reaches)
geoms_buff <- pbapply::pblapply(1:nrow(reaches), \(i) {
    # maximum width is half the reach length. 
    width = min(as.numeric(st_length(reaches[i, ])/2), reaches$width[i]+2*sqrt(reaches$width_var[i]))
    st_buffer(geoms[i], width)
}) 
# combine into one sfc object
geoms_buff <- Reduce(c, geoms_buff)
reaches$buff_geom <- geoms_buff

# Create a polygon within each buffer which includes only 
# water pixels from GRWL water mask which are >20m from shore 
# (to help with edge effects from atmospheric correction)


get_polygon <- function(i) {
    reach_id = reaches[i,]$reach_id
    r <- rast(vect(geoms_buff[i]), res = 10)
    r <- rasterize(vect(geoms_buff[i]), r)
    # Make distance rasters for all overlapping polygons
    neighbors <- st_intersection(reaches, geoms_buff[i, ] )
    for (g in st_geometry(neighbors)) {
        d <- rasterize(vect(g), r)
        r <- c(r, distance(d))
    }
    # Make raster with index of closest reach
    r <- r[[-1]]
    # way faster than apply(x, 1, which.min)
    closest_vect <- max.col(-as.matrix(r))
    closest_r <- rast(d)
    values(closest_r) <- neighbors$reach_id[closest_vect]
    rm(closest_vect)
    # Make mask of water pixels closest to reach
    water_mask <- crop(grwl_mask, project(ext(r), closest_r, grwl_mask))
    water_mask <- project(water_mask, closest_r, method = "near")
    # 256 = No Data
    # 255 = River
    # 180 = Lake/reservoir
    # 126 = Tidal rivers/delta
    # 86  = Canal
    reach_mask <- (water_mask %in% c(86, 180, 255, 126)) & (closest_r == reach_id)
    # make polygon
    p <- as.polygons(reach_mask)
    p_ <- p[p[[1]][[1]] == 1, ] #Only include "true" polygon (it also makes an inverse of "false" values)
}

polygon_list <- pbapply::pblapply(1:nrow(reaches), get_polygon)
polygon_list <- pbapply::pblapply(polygon_list, st_as_sf)
water_polygons <- Reduce(rbind, polygon_list)
write_sf(water_polygons, file.path(basin_dir, "water_polygons.gpkg"))

# ## make a quick figure for general exam thing
# # first need a bbox to include
# # 46.520432, -119.450801, 45.896930, -118.718603
# bounds <- data.frame(lon = c(-118.718603, -119.450801),
#                    lat = c(46.520432, 45.896930)) |>
#   st_as_sf(coords = c("lon", "lat"), 
#            crs = 4326) |>
#   st_bbox() |>
#   st_as_sfc()
# d <- st_intersection(water_polygons, st_transform(bounds, st_crs(water_polygons)))
# ggplot(d) + 
#     geom_sf(fill = rep_len(c("green", "blue"), 25))
