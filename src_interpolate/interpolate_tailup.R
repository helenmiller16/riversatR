library(data.table)
library(Matrix)
library(RTMB)
source(here::here("src_interpolate/functions.R"))
reaches <- sf::read_sf("data/external/mississippi_small/reaches.gpkg")
reaches <-dplyr::rename(reaches, 
                        reach_id_up = "rch_id_up", 
                        reach_id_down = "rch_id_dn")
reaches$reach_len <- reaches$reach_len / 1000 # convert to km
reflectance_dt <- fread("data/external/mississippi_small/landsat_2015_10_01-2016_10-01.csv")
reflectance_dt[, ndti := (red - green)/(red + green)]
variables = c("ndti")
graph <- make_graph(reaches)
data <- make_data(reflectance_dt, 
                  graph, 
                  variables = variables,
                  time_interval = 8, 
                  weight_variable = "facc")
data$y_n_t_v <- data$y_n_t_v[, 22:24, ,drop = FALSE]
prod(dim(data$y_n_t_v)) # total number of variables to estimate. 
params <- list(log_theta = -5, 
               log_gamma1 = 0, 
               log_gamma2 = 0, 
               log_sigma = 0, 
               mu = 0, 
               z_n_t = rep(0, prod(dim(data$y_n_t_v))))

f <- objective_function(tailup_no_covariates, 
                        data, 
                        params)

obj <- MakeADFun(f, params, random = "z_n_t")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)

rand <- summary(sdr, "random")
fixed <- summary(sdr, "fixed")
yhat <- rand[, "Estimate"] + fixed["mu", "Estimate"]
plot(data$y_n_t_v[, 1, 1], yhat)
reaches$obs <- data$y_n_t_v[, 1, 1]
reaches$yhat <- yhat

library(ggplot2)
ggplot(reaches) + 
  geom_sf(aes(color = obs)) +
  scale_color_continuous(limits = c(-.5, .5), 
                         palette = "viridis", 
                         na.value = "transparent") + 
  theme_minimal()

ggplot(reaches) + 
  geom_sf(aes(color = yhat)) +
  scale_color_continuous(limits = c(-.5, .5), 
                         palette = "viridis", 
                         na.value = "transparent") + 
  theme_minimal()

# 
# # plot estimates
# z_t <-  data.table(x = rand[, "Estimate"])
# z_t$time <- days
# z_t$se <- rand[, 'Std. Error']



# FOR TESTING
if (FALSE) {
  library(sf)
  library(data.table)
  library(terra)
  reaches <- sf::read_sf("data/external/mississippi_small/reaches.gpkg")
  reaches <-dplyr::rename(reaches, 
                          reach_id_up = "rch_id_up", 
                          reach_id_down = "rch_id_dn")
  reflectance_dt <- fread("data/external/mississippi_small/landsat_2015_10_01-2016_10-01.csv")
  reaches$reach_len <- reaches$reach_len / 1000 # convert to km
  variables = c("ndti")
  reflectance_dt[, ndti := (red-green)/(red+green)]
  graph <- make_graph(reaches)
  data <- make_data(reflectance_dt, 
                    graph, 
                    variables = variables,
                    time_interval = 8, 
                    weight_variable = "facc")
  # 
  nodes_dt <- as.data.table(sfnetworks::activate(graph, "nodes"))
  # should be in order already 
  reaches$node_id <- nodes_dt[.(reaches$reach_id), node_id, on = "reach_id"]
  reaches <- reaches[order(reaches$node_id), ]
  obs <- apply(data$y_n_t_v, 1, \(x) sum(!is.na(x))/length(variables))
  reaches$obs <- obs
  reaches <- terra::vect(reaches)
  plot(reaches, 'obs')
  # PARAMS? 
  
  data$y_n_t_v <- data$y_n_t_v[, 22, ,drop = FALSE]
  prod(dim(data$y_n_t_v)) # total number of variables to estimate. 
  params <- list(log_theta = -6, 
                 log_gamma1 = 0, 
                 log_gamma2 = 0, 
                 log_sigma = -5, 
                 mu = 0, 
                 z_n_t = rep(0, prod(dim(data$y_n_t_v))))
  Q <- make_prec(data, params, tailup = TRUE)
  C <- solve(Q)
  G <- make_prec(data, params, tailup = TRUE, precision = F)
  image(G[1:30, 1:30])
  
  sim <- rmvnorm_prec(rep(0, ncol(Q)), prec = Q)
  reaches$sim1 <- sim
  plot(reaches, 'sim1', breaks = 20)
  # idx <- rep(1:3, each = dim(data$y_n_t_v)[1])
  # idx <- rep(1:3, length.out= nrow(Q))
  reaches$obs1 <- data$y_n_t_v[, 1, 1]
  
  plot(reaches, 'obs1', range = c(-.1, .1), breaks = 20)
}