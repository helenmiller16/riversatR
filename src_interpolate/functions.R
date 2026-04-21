# This may be helpful: 
# https://github.com/James-Thorson-NOAA/graphical_mixed_model/blob/main/trait_functions.R

require(data.table)
require(Matrix)
require(RTMB)


# function from textbook
rmvnorm_prec <-
  function( mu, # estimated fixed and random effects
            prec, # estimated joint precision
            n.sims = 1) {
    
    require(Matrix)
    # Simulate values
    z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    # Q = t(P) * L * t(L) * P
    L = Cholesky(prec, super=TRUE)
    # Calcualte t(P) * solve(t(L)) * z0 in two steps
    z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
    z = solve(L, z, system = "Pt") # z = Pt    * z
    return(mu + as.matrix(z))
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
  
  nodes <- sf::st_as_sf(sfnetworks::activate(graph, "nodes"))
  nodes_dt <- data.table::as.data.table(nodes)
  edges_dt <- data.table::as.data.table(sfnetworks::activate(graph, "edges"))
  ref$node_id <- nodes_dt[.(ref$reach_id), node_id, on = "reach_id"]
  n <- nrow(nodes)
  v <- length(variables)
  t <- trunc(as.numeric(max(ref$date) - min(ref$date))/time_interval) +1
  dims <- c(n, t, v)
  
  # get time step for each obs. 
  # first `time_interval` days go in time_step 1 etc. 
  date <- seq(min(ref$date), max(ref$date), by = 1)
  interval_date <- seq(min(ref$date), max(ref$date), by = time_interval)
  time_step <- rep(1:length(interval_date), each = time_interval, length.out = length(date))
  time_step <- data.table(date, time_step)
  ref <- merge(ref, time_step)
  # ref[, date_diff := time_step[.(ref$date), on = "date"]$date - ref$date]
  ref[, date_diff := as.numeric(date - interval_date[time_step])]
  ref_time_step <- ref[, lapply(.SD, mean, na.rm = TRUE), .(time_step, reach_id)]
  ref_date_time_step <- ref[, .(date_diff = unique(date_diff)), .(reach_id)]
  # 
  # data_list <- lapply(dates, \(d) {
  #   days <- seq(d, by = 1, length.out = time_interval)
  #   ref[date %in% days][
  #     , lapply(.SD, mean, na.rm = TRUE), node_id][
  #       .(1:n), on = "node_id"][,..variables]
  # })
  
  data_list <- lapply(variables, \(v) {
    d <- dcast(ref_time_step, node_id ~ time_step, 
               value.var = v,)
    d <- d[.(1:n), on = "node_id"]
    d$node_id <- NULL
    as.matrix(d)
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
    flow <- nodes_dt[[weight_variable]]
  }
  
  weights_e <- sapply(1:nrow(edges_dt), \(e) {
    from_node <- edges_dt$from[e]
    to_node <- edges_dt$to[e]
    flow[from_node]/sum(flow[edges_dt[to == to_node]$from])
  })
  
  list(
    y_n_t_v = array(data = unlist(data_list), dim = dims),
    from_e = edges_dt$from,
    to_e = edges_dt$to,
    dist_e = edges_dt$length,
    weights_e = weights_e,
    time_step_date_map = ref_date_time_step
    # weights_mx = weights_mx
  )
}


make_prec_tailup <- function(data, 
                      log_gamma1 = 0, 
                      log_gamma2 = 0, 
                      log_gamma3 = 0, 
                      log_theta  = 0,
                      type = c("space", "spacetime", "time")) {
  
  N <- dim(data$y_n_t_v)[1]
  E <- length(data$from_e)
  to_e <- data$to_e
  from_e <- data$from_e
  dist_e <- data$dist_e
  theta <- exp(log_theta)
  sources <- which(!1:N %in% to_e)
  
  path_e <- exp(-theta*dist_e) * data$weights_e
  
  # Path_tailup <- AD(Matrix(data = 0, nrow = N, ncol = N, sparse = TRUE), force = TRUE)
  Path_tailup <- AD(sparseMatrix(i = to_e, 
                                 j = from_e, 
                                 x = 0, 
                                 dims = c(N,N)), force = TRUE)
  Path_tailup[cbind(to_e, from_e)] = path_e
  v_tailup <- AD(rep(0, N))
  
  for (i in 1:E) {
    # variance: weighted average for each x
    w <- 1- (exp(-2*theta*dist_e[i]) * data$weights_e[i])
    v_tailup[to_e[i]] <- v_tailup[to_e[i]] + w
  }
  # apply(var_mx, 1, sum)
  # set variance for initial conditions
  for (s in sources) {
    v_tailup[s] <- 1
  }
  t <- dim(data$y_n_t_v)[2]
  I_t <- AD(Diagonal(t))
  v_t <- rep(1, t)
  n_pjoint <- N * t
  
  Lag1 <- AD(Matrix(data = 0, nrow = t, ncol = t, sparse = TRUE))
  for (i in 1:(t-1)) {
    Lag1[1+i, i] <- 1
  }
  
  if ("spacetime" %in% type) {
    # dgmrf fails with more than 7 time steps if we try to 
    # evaluate this gamma2 term for some reason...? 
    gamma1 <- AD(Diagonal(n_pjoint), force = TRUE)
    gamma1[cbind(1:n_pjoint, 1:n_pjoint)] <- exp(log_gamma1)
    # gamma3 <- AD(Diagonal(n_pjoint), force = TRUE)
    # gamma3[cbind(1:n_pjoint, 1:n_pjoint)] <- exp(log_gamma3)
    
    I_ptailup <- AD(Diagonal(N), force = TRUE)
    
    P_spacetime <- gamma1 %*% kronecker(Path_tailup, Lag1)# spatiotemporal
  } else {
    P_spacetime <- AD(Matrix(data = 0, nrow= n_pjoint, ncol = n_pjoint), force = TRUE)
  }
  
  if ("space" %in% type) {
    gamma2 <- AD(Diagonal(n_pjoint), force = TRUE)
    gamma2[cbind(1:n_pjoint, 1:n_pjoint)] <- exp(log_gamma2)
    P_space <- gamma2 %*% kronecker(Path_tailup, I_t)
  } else {
    P_space <- AD(Matrix(data = 0, nrow= n_pjoint, ncol = n_pjoint), force = TRUE)
  }
  
  if ("time" %in% type) {
    I_n <- AD(Diagonal(N), force = TRUE)
    gamma3 <- AD(Diagonal(n_pjoint), force = TRUE)
    gamma3[cbind(1:n_pjoint, 1:n_pjoint)] <- exp(log_gamma3)
    P_time <- gamma3 %*% kronecker(I_n, Lag1)
  } else {
    P_time <- AD(Matrix(data = 0, nrow= n_pjoint, ncol = n_pjoint), force = TRUE)
  }
  
  P_joint = P_spacetime + P_space + P_time
  
  v <- kronecker(v_tailup, v_t)
  v_inv <- 1/v
  V_inv <- AD(Diagonal(n_pjoint), force = TRUE)
  V_inv[cbind(1:n_pjoint, 1:n_pjoint)] <- v_inv
  
  I_pjoint <- AD(Diagonal(n_pjoint), force = TRUE)
  t(I_pjoint-P_joint) %*% V_inv %*% (I_pjoint-P_joint)
}


tailup_no_covariates <- function(data, params) {
  # data = list with : 
  # y_n_t_v: array where dim 1 is location, dim 2 is time, dim 3 is variable.
  # from_e length e (edges)
  # to_e length e (edges)
  # dist_e length e (edges)
  # weights_mx 
  # params = list with 
  # log_theta
  # log_gamma1
  # log_gamma2
  # log_sigma
  # mu
  # z_n_t_f # matrix
  data$y_n_t_v <- OBS(data$y_n_t_v)
  jnll_random <- 0
  jnll_fixed <- 0
  
  Q <- make_prec_tailup(data, 
                        log_gamma1 = params$log_gamma1, 
                        log_gamma2 = params$log_gamma2,
                        log_theta = params$log_theta, 
                        type = c("space", "spacetime"))
  jnll_random <-  - dgmrf(params$z_n_t, mu = 0, Q = Q, log = TRUE)
  
  x_n_t <- params$mu + params$z_n_t
  jnll_fixed <- jnll_fixed - sum(dnorm(as.vector(data$y_n_t_v), 
                                       mean = x_n_t, 
                                       sd = exp(data$log_sigma), 
                                       log = TRUE), na.rm = TRUE)
  

  REPORT(jnll_random)
  REPORT(jnll_fixed)
  REPORT(Q)
  REPORT(x_n_t)
  ADREPORT(x_n_t)
  jnll_random + jnll_fixed
}

tailup_mean <- function(data, params) {
  # data = list with : 
  # y_n_t_v: array where dim 1 is location, dim 2 is time, dim 3 is variable.
  # from_e length e (edges)
  # to_e length e (edges)
  # dist_e length e (edges)
  # weights_mx 
  # params = list with 
  # log_theta
  # log_gamma1
  # log_gamma2
  # log_sigma
  # mu
  # w_n
  # z_n_t # 
  data$y_n_t_v <- OBS(data$y_n_t_v)
  jnll_random <- 0
  jnll_fixed <- 0
  t <- dim(data$y_n_t_v)[2]
  N <- dim(data$y_n_t_v)[1]
  
  
  Q_space <- make_prec_tailup(data, 
                              log_theta = params$log_theta,
                              type = "space")[1:N, 1:N]
  jnll_random <- jnll_random - dgmrf(params$w_n, 
                                     mu = 0,
                                     Q = Q_space, 
                                     log = TRUE, 
                                     scale = 1/exp(params$log_gamma2))
 
  
  
  x_n_t <- params$mu + rep(params$w_n, t) 
  jnll_fixed <- jnll_fixed - sum(dnorm(as.vector(data$y_n_t_v), 
                                       mean = x_n_t, 
                                       sd = exp(params$log_sigma), 
                                       log = TRUE), na.rm = TRUE)
  
  
  REPORT(jnll_random)
  REPORT(jnll_fixed)
  REPORT(Q_space)
  REPORT(x_n_t)
  REPORT(params$w_n)
  ADREPORT(x_n_t)
  jnll_random + jnll_fixed
}

tailup_iter_time <- function(data, params) {
  # data = list with : 
  # y_n_t_v: array where dim 1 is location, dim 2 is time, dim 3 is variable.
  # from_e length e (edges)
  # to_e length e (edges)
  # dist_e length e (edges)
  # weights_mx 
  # params = list with 
  # log_theta
  # log_gamma1
  # log_gamma2
  # log_sigma
  # log_rhoT
  # mu
  # w_n
  # z_n_t # 
  data$y_n_t_v <- OBS(data$y_n_t_v)
  jnll_random <- 0
  jnll_fixed <- 0
  t <- dim(data$y_n_t_v)[2]
  N <- dim(data$y_n_t_v)[1]
  z_n_t <- array(params$z_n_t, dim = dim(data$y_n_t_v))
  w_n <- params$w_n
  
  Q_space <- make_prec_tailup(data, 
                              log_theta = params$log_theta,
                              type = "space")[1:N, 1:N]
  jnll_random <- jnll_random - dgmrf(w_n,
                                     mu = 0,
                                     Q = Q_space,
                                     scale = exp(params$log_gamma2),
                                     log = TRUE)
  # Q_spacetime <- make_prec_tailup(data, params, type = "spacetime")
  # jnll_random <- jnll_random - dgmrf(params$z_n_t, 
  #                                    mu = 0, 
  #                                    Q = Q_spacetime, 
  #                                    log = TRUE)
  
  for (i in 1:t) {
    if (i==1) {
      jnll_random <- jnll_random - dgmrf(z_n_t[, i], 
                                         mu = 0, 
                                         Q=Q_space,
                                         scale = (exp(params$log_gamma1)),
                                         log = TRUE)
    } else {
      jnll_random <- jnll_random - dgmrf(z_n_t[, i] - exp(params$log_rhoT) * z_n_t[, i-1], 
                                         mu = 0, 
                                         Q=Q_space,
                                         scale = exp(params$log_gamma1),
                                         log = TRUE)
    }
  }
  
  
  x_n_t <- params$mu + rep(w_n, t) + params$z_n_t
  jnll_fixed <- jnll_fixed - sum(dnorm(as.vector(data$y_n_t_v),
                                       mean = x_n_t,
                                       sd = exp(params$log_sigma),
  log = TRUE), na.rm = TRUE)
  
  # jnll_fixed <- -sqrt(sum((data$y_n_v_t - x_n_t)^2))
  
  
  REPORT(jnll_random)
  REPORT(jnll_fixed)
  REPORT(Q_space)
  REPORT(x_n_t)
  REPORT(z_n_t)
  REPORT(w_n)
  ADREPORT(x_n_t)
  jnll_random + jnll_fixed
}

objective_function <- function(fun, data, params) {
  function(params) {
    fun(data, params)
  }
}