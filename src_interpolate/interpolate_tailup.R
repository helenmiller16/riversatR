# FOR TESTING
if (FALSE) {
  reaches <- read_sf("data/external/columbia/reaches.gpkg")
  reaches <-dplyr::rename(reaches, 
                          reach_id_up = "rch_id_up", 
                          reach_id_down = "rch_id_dn")
}


# Reaches is a sf object with a row for each 
# river reach, and the columns
# reach_id
# reach_id_up
# reach_id_down
# if more than 1 upstream or downstream reach then 
# from_id or to_id may have multiple ids separated with a space
# (This is the format of the SWOT dataset)
make_graph <- function(reaches) {
  
  edges <- as.data.table(reaches)
  
  # Create "nodes" table by looping through all the 
  # edges and adding nodes as needed. add node geometry
  # from the ends of the edges. It won't be perfect 
  # (e.g. edges may not necessarily line up with nodes) 
  # but topology should be correct. 
  n <- 0
  edges$from <- 0
  edges$to <- 0
  edge_geometry <- lapply(edges$geom, \(g) g[[1]])
  node_geometry <- list()
  for (i in 1:nrow(edges)) {
    edge <- edges[i, ]
    from_e <- as.numeric(strsplit(edge$reach_id_up, " ")[[1]])
    to_e <- as.numeric(strsplit(edge$reach_id_down, " ")[[1]])
    
    if (edge$from == 0) {
      if (any(edges[reach_id %in% from_e]$to_node>0)){
        if (length(from_e) > 1) {
          stop("length(from_e) > 1")
        }
        edges[i]$from <- edges[reach_id %in% from_e]$to_node
      } else {
        n <- n + 1
        edges[i]$from <- n
        edges[reach_id %in% from_e]$to <- n
        node_geometry <- c(node_geometry, list(st_point(edge_geometry[[i]][1, ])))
      }
    }
    
    if (edge$to == 0) {
      if (sum(edges$reach_id %in% to_e) > 0) {
        if (any(edges[reach_id %in% to_e]$from>0)) {
          if (length(to_e) > 1) {
            stop("length(to_e) > 1")
          }
          edges[i]$to <- edges[reach_id %in% to_e]$from
        } else {
          n <- n + 1
          edges[i]$to <- n
          edges[reach_id %in% to_e]$from <- n
          node_geometry <- c(node_geometry, list(st_point(edge_geometry[[i]][nrow(edge_geometry[[i]]), ])))
        }
      } else {
        n <- n + 1
        edges[i]$to <- n
        node_geometry <- c(node_geometry, list(st_point(edge_geometry[[i]][nrow(edge_geometry[[i]]), ])))
      }
    }
  }
  
  # Make nodes shapefile
  N <- max(c(edges$from, edges$to))
  nodes <- st_sf(node_id = 1:N, geometry = st_sfc(node_geometry))
  
  # Check for any edges which don't line up with nodes. 
  # First check for nodes which are far apart
  st_crs(nodes) <- st_crs(reaches)
  nodes_proj <- st_transform(nodes, "EPSG:5070")
  edges_distance <- st_distance(nodes[edges$from,], nodes[edges$to,], by_element = TRUE)
  edges$geom <- st_sfc(lapply(edges$geom, st_cast, to = "LINESTRING"))
  edges <- st_as_sf(edges, crs = st_crs(reaches))
  sfnetworks::sfnetwork(nodes, 
                                   edges, 
                                   directed = TRUE, 
                                   force = TRUE)

}

# FOR TESTING
if (FALSE) {
  include <- st_read("data/external/columbia/region_include.gpkg")
  reaches <- read_sf("data/external/columbia/reaches.gpkg")
  reaches <-dplyr::rename(reaches, 
                          reach_id_up = "rch_id_up", 
                          reach_id_down = "rch_id_dn")
  include <- st_transform(include, st_crs(reaches))
  reaches <- reaches[include, ]
  ref <- fread("data/external/conus/reflectance.tsv")
  ref <- ref[reach_id %in% reaches$reach_id]
  reaches_ <- reaches[reaches$reach_id %in% ref$reach_id, ]
  
  reflectance_dt <- ref
  variables = c("coastal", 
                "blue", 
                "green", 
                "red", 
                "red_edge1", 
                "red_edge2", 
                "red_edge3", 
                "red_edge4", 
                "nir")
  
  
}

make_data <- function(
  # data.table with columns: 
  # reach_id lining up with reaches$reach_id
  # date (POSIXct or Datetime object, coverted to Datetime)
  # variables listed in 'variables'
  reflectance_dt, 
  graph, 
  variables = c("coastal", 
                "blue", 
                "green", 
                "red", 
                "red_edge1", 
                "red_edge2", 
                "red_edge3", 
                "red_edge4", 
                "nir")) {
  # We're going to start by trying daily predictions. 
  ref <- copy(reflectance_dt)
  ref[, date := as.Date(date)]
  
  nodes <- st_as_sf(sfnetworks::activate(graph, "nodes"))
  edges <- st_as_sf(sfnetworks::activate(graph, "edges"))
  edges_dt <- as.data.table(edges)
  setindex(edges_dt, "reach_id")
  # Give each node the value of the downstream edge
  ref$node <- edges_dt[.(ref$reach_id), from, on = "reach_id"]
  
  n <- nrow(nodes)
  v <- length(variables)
  t <- as.numeric(max(ref$date) - min(ref$date)) +1
  dims <- c(n, v, t)
  
  dates <- seq(min(ref$date), max(ref$date), by = "day")
  
  data_list <- lapply(dates, \(d) {
    ref[date == d][.(1:n), on = "node"][, ..variables]
  })
  
  # NEED TO GET FLOW STILL
  
  array(data = unlist(data_list))
}


tailup_no_covariates <- function(data, params) {
  # data = list with : 
    # y_n_v_t: array where dim 1 is location, dim 2 is variable, dim 3 is time.
    # from_e length n-1
    # to_e length n-1
    # dist_e length n-1
    # flow_n length n
  # params = list with 
  
}

objective_function <- function(fun, data, params) {
  function(data) {
    fun(data, params)
  }
}