colors <- function(x, palette = "viridis", reverse = FALSE, n = 10, log=  FALSE, return_breaks = FALSE) {
  
  if (palette == "viridis") {
    cols <- c("#440154", "#31688e","#35b779","#fde725")
  } else if (palette == "inferno") {
    cols <- c("#fcffa4", "#ed6925", "#781c6d", "#000004")
  } else if (palette == "brown") {
    cols <- c("#E7C7A8FF", "#874D4AFF", "#3B141AFF", "#27170BFF")
  } else if (palette == "green") {
    cols <- paletteer::paletteer_d("ggsci::green_material")
  } else if (palette == "yellow_green") {
    cols <- paletteer::paletteer_d("MoMAColors::Alkalay2")
  } else if (palette == "teal") {
    cols <- c("#B2DFDBFF", "#00897BFF", "#002620")
  } else if (palette == "teal_colorful") {
    cols <- c("#EDEF5DFF", "#39AB7EFF",  "#012423",)
  }
  
  if (reverse) cols <- rev(cols)
  
  pal <- colorRampPalette(cols)(n)
  if (log) {
    x_ <- log(x)
    x.inf <- !is.na(x_) & x_ == -Inf
    x_[x.inf] <- NA
    x_[x.inf] <- min(x_, na.rm = TRUE)
  } else x_ <- x
  
  
  breaks <- seq(min(x_, na.rm = TRUE), max(x_,na.rm = TRUE), length = n+1)
  if (return_breaks) {
    if (log) breaks <- exp(breaks)
    return(list(breaks = breaks, 
                cols = pal))
  }
  
  b <- .bincode(x_, breaks = breaks, include.lowest = TRUE)
  pal[b]
  
}


color_legend <- function(x, x_ = "topleft", y = NULL, title = "",  palette = "viridis", reverse = FALSE, n = 10, log=  FALSE) {
  p <- colors(x, palette, reverse, n, log, return_breaks = TRUE)
  labels <- c(format(p$breaks[1], digits = 2), rep("", n-2), format(p$breaks[n+1], digits = 2))
  legend(x_, 
         legend = format(labels, digits = 2), 
         fill = p$cols, 
         title = title, 
         bty = "n", 
         y.intersp = .5)
}


# For fetching PC Sentinel data ----- 

# Get a list of rasts
get_rast <- function(date, band, items) {
  
  if (!is.null(date)) {
    features <- items$features[items_datetime(items) == date]
  } else {
    features <- items$features
  }
  
  
  rast_list <- lapply(features, \(f) rast(paste0("/vsicurl/", f$assets[[band]]$href)))
  # SpatRasterCollection
  rast_sprc <- sprc(rast_list)
  # algo = 3 creates a virtual raster
  r <- merge(rast_sprc, algo = 3)
  if (!is.null(date)) {
    time(r) <- as.POSIXct(date)
  }
  r
}

# return data.table with a row for each non-NA value
extract_data <- function(rast, # rast with 'time' attribute and one layer corresponding to 'band'
                         locations, # geom of locations to extract for
                         dt_id, # data.table with identifying info for locations
                         band # band name
) {
  dt <- copy(dt_id)
  dt[['date']] <- time(rast)
  dt[['band']] <- band
  date <- time(rast)
  # handle data capture error
  ex <- try(extract(rast, locations)[, 2])
  if ("try_error" %in% class(ex)) {
    warning("Failed to fetch data for ", names(rast))
    dt_id[['value']] <- NA
    return(dt_id[0, ])
  }
  dt[['value']] <- ex
  dt[ !is.na(dt[['value']]) ]
}



get_data <- function(band, items) {
  
  
  # separate out by crs
  items_df <- items_as_tibble(items)
  
  get_data_list <- function(crs) {
    
    i <- items |> 
      items_filter(properties$`proj:epsg` == crs)
    
    i_df <- items_as_tibble(i)
    
    rast_list <- lapply(unique(i_df$datetime), 
                        get_rast, 
                        band = band, 
                        items = i)
    
    data_list <- lapply(rast_list, extract_data, 
                        locations =  st_transform(soc_data, paste0("EPSG:", crs)), 
                        dt_id = dt_id, 
                        band = band)
    
    rbindlist(data_list)
  }
  
  data_list <- lapply(unique(items_df$`proj:epsg`), get_data_list) 
  
  rbindlist(data_list)
}