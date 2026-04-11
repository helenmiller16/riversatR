require(data.table)
# Reaches is a sf object with a row for each 
# river reach, and the columns
# reach_id
# reach_id_up
# reach_id_down
# if more than 1 upstream or downstream reach then 
# from_id or to_id may have multiple ids separated with a space
# (This is the format of the SWOT dataset)
make_graph <- function(reaches) {
  nodes <- data.table::as.data.table(reaches)
  nodes$node_id <- 1:nrow(nodes)
  nodes$geom <- sf::st_line_interpolate(sf::st_transform(nodes$geom, "EPSG:5070"), .5, TRUE)
  
  # Make edges
  from <- lapply(nodes$reach_id_up, \(e) as.numeric(strsplit(e, ' ')[[1]]))
  to   <- nodes$reach_id
  
  from <- c(from, nodes$reach_id)
  to   <- c(to, lapply(nodes$reach_id_down, \(e) as.numeric(strsplit(e, ' ')[[1]])))
  
  edge_list <- mapply(\(f, t) {
    data.table::data.table(expand.grid(f, t))
  }, from, to, SIMPLIFY = FALSE)
  edges <- data.table::rbindlist(edge_list)
  names(edges) <- c("from", "to")
  edges <- edges[!duplicated(edges)]
  
  # remove edges which come from reaches not in the network
  edges <- edges[from %in% nodes$reach_id]
  edges <- edges[to %in% nodes$reach_id]
  
  # Calculate length of each edge. Since the go from the center of two
  # reaches, it is the mean length of those two reaches
  edges$from_len <- nodes[.(edges$from), on = "reach_id"]$reach_len
  edges$to_len <- nodes[.(edges$to), on = "reach_id"]$reach_len
  edges[, length := (from_len + to_len)/2]
  edges$from <- nodes[.(edges$from), on = "reach_id"]$node_id
  edges$to <- nodes[.(edges$to), on = "reach_id"]$node_id
  
  sfnetworks::sfnetwork(
    sf::st_as_sf(nodes), 
    edges, 
    directed = TRUE)
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
                "nir08", 
                "swir16", 
                "swir22"), 
  # Time step/interval. If there are multiple observations at a site
  # within a time interval take the mean. 
  time_interval = 8, 
  # weight_variable is a column in the graph nodes object to be used for 
  # weighting inputs at confluences (e.g. contributing area). 
  # If NA, weights will be the number of upstream confluences. 
  # TODO: allow this to vary over time (and be related to discharge)
  weight_variable = "facc" # facc is flow accumulation in SWOT v17b
  ) {
  ref <- copy(reflectance_dt)
  ref[, date := as.Date(date)]
  
  nodes <- st_as_sf(sfnetworks::activate(graph, "nodes"))
  nodes_dt <- data.table::as.data.table(nodes)
  edges_dt <- data.table::as.data.table(sfnetworks::activate(graph, "edges"))
  ref$node_id <- nodes_dt[.(ref$reach_id), node_id, on = "reach_id"]
  n <- nrow(nodes)
  v <- length(variables)
  t <- trunc(as.numeric(max(ref$date) - min(ref$date))/time_interval) +1
  dims <- c(n, v, t)
  
  dates <- seq(min(ref$date), max(ref$date), by = time_interval)
  
  data_list <- lapply(dates, \(d) {
    days <- seq(d, by = 1, length.out = time_interval)
    ref[date %in% days][
      , lapply(.SD, mean, na.rm = TRUE), node_id][
        .(1:n), on = "node_id"][,..variables]
  })
  
  
  
  if (is.na(weight_variable)) {
    # get flow and sources
    sources <- which(!1:nrow(nodes) %in% edges_dt$to)
    # use flow accumulation (facc) 
    source_flow <- sapply(sources, \(s) nodes_dt[node_id == s]$facc)
    flow <- rep(0, n)
    flow[sources] <- source_flow
    # This step takes a while for large graphs
    for (source in sources) {
      # get list of all downstream points
      connected <- sapply(1:nrow(nodes), 
                          \(x) ifelse(source == x, 0, igraph::edge_connectivity(graph, source, x)))
      # add source flow to all downstream points
      flow <- flow + unlist(connected)*flow[source]
    }
  } else {
    flow <- nodes_dt[, ..weight_variable]
  }
  
  
  list(
    y_n_v_t = array(data = unlist(data_list), dim = dims),
    from_e = edges_dt$from,
    to_e = edges_dt$to,
    dist_e = edges_dt$length,
    flow_n = flow
  )
}


# FOR TESTING
if (FALSE) {
  library(sf)
  library(data.table)
  reaches <- sf::read_sf("data/external/mississippi_small/reaches.gpkg")
  reaches <-dplyr::rename(reaches, 
                          reach_id_up = "rch_id_up", 
                          reach_id_down = "rch_id_dn")
  reflectance_dt <- fread("data/external/mississippi_small/landsat_2015_10_01-2016_10-01.csv")
  variables = c("coastal", 
                "blue", 
                "green", 
                "red", 
                "nir08", 
                "swir16", 
                "swir22")
  graph <- make_graph(reaches)
  data <- make_data(reflectance_dt, 
                    graph, 
                    variables = variables,
                    time_interval = 8, 
                    weight_variable = "facc")
  # PARAMS? 
  
}


tailup_no_covariates <- function(data, params) {
  # data = list with : 
    # y_n_v_t: array where dim 1 is location, dim 2 is variable, dim 3 is time.
    # from_e length e (edges)
    # to_e length e (edges)
    # dist_e length e (edges)
    # flow_n length n (nodes)
  # params = list with 
  
}


objective_function <- function(fun, data, params) {
  function(data) {
    fun(data, params)
  }
}

main <- function(reachs_sf, reflectance_dt) {
  graph <- make_graph(reaches_sf)
  data <- make_data(reflectance_dt, graph)
}