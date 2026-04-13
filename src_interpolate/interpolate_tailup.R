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
y_n_t_v <- data$y_n_t_v
# y_n_t_v[y_n_t_v < -1] <- NA
# y_n_t_v[y_n_t_v > 1] <- NA
n_fits <- ceiling(dim(y_n_t_v)[2]/4)
starts <- seq(1, dim(y_n_t_v)[2], by = 3)
for (start in starts) {
  end <- start + 5
  data$y_n_t_v <- y_n_t_v[, 1:5, ,drop = FALSE]
  params <- list(log_theta = -5, 
                 log_gamma1 = 0,
                 log_gamma2 = 0,
                 log_gamma3 = 0,
                 log_sigma = 0, 
                 mu = 0, 
                 z_n_t = rnorm(prod(dim(data$y_n_t_v))),
                 w_n = rnorm(dim(data$y_n_t_v)[1]))
  
  f <- objective_function(tailup_no_covariates, 
                          data, 
                          params)
  # f2 <- objective_function(tailup_iter_time, 
  #                          data, 
  #                          params)
  # params2 <- copy(params)
  # params2$z_n_t <- params$z_n_t-rep(params$w_n, dim(data$y_n_t_v)[2])

  obj <- MakeADFun(f, params, random = "z_n_t")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  saveRDS(opt, paste0("data/external/mississippi_small/results/opt_", start, ".rds"))
  r <- obj$report()
  yhat <- r$x_n_t
  yhat <- array(yhat, dim = dim(data$y_n_t_v))
  saveRDS(yhat, paste0("data/external/mississippi_small/results/yhat_", start, ".rds"))
  sdr <- sdreport(obj)
  saveRDS(summary(sdr, "random"), paste0("data/external/mississippi_small/results/sdr_random_", start, ".rds"))
  saveRDS(summary(sdr, "fixed"), paste0("data/external/mississippi_small/results/sdr_fixed_", start, ".rds"))
  rm(sdr)
  rm(obj)
  rm(opt)
}
# data$y_n_t_v <- y_n_t_v[, 1:5, ,drop = FALSE]
# params <- list(log_theta = -5, 
#                log_gamma1 = 0,
#                log_gamma2 = 0,
#                log_gamma3 = 0,
#                log_sigma = 0, 
#                mu = 0, 
#                z_n_t = rnorm(prod(dim(data$y_n_t_v))))
# 
# f <- objective_function(tailup_no_covariates, 
#                         data, 
#                         params)
# f(params)
# obj <- MakeADFun(f, params, random = "z_n_t")
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# saveRDS(opt, "data/external/mississippi_small/results/opt.rds")
# r <- obj$report()
# yhat <- r$x_n_t
# yhat <- array(yhat, dim = dim(data$y_n_t_v))
# saveRDS(yhat, "data/external/mississippi_small/results/yhat.rds")
# plot(data$y_n_t_v[, , 1], yhat)
# abline(0,1, col = "blue")
# reaches$obs1 <- data$y_n_t_v[, 1, 1]
# reaches$obs2 <- data$y_n_t_v[, 2, 1]
# reaches$obs3 <- data$y_n_t_v[, 3, 1]
# reaches$yhat1 <- yhat[, 1, 1]
# reaches$yhat2 <- yhat[, 2, 1]
# reaches$yhat3 <- yhat[, 3, 1]
# 
# library(ggplot2)
# library(patchwork)
# p <- lapply(c("obs1", "obs2", "obs3", 
#               "yhat1", "yhat2", "yhat3"), \(x) {
#                 ggplot(reaches) + 
#                   geom_sf(aes(color = .data[[x]])) +
#                   scale_color_continuous(limits = c(-.5, .5), 
#                                          palette = "viridis", 
#                                          na.value = "transparent") + 
#                   theme_minimal()
#               })
# 
# 
# fig <- (p[[1]] + p[[2]] + p[[3]]) / (p[[4]] + p[[5]] + p[[6]])
# pdf("fitted.pdf", width = 11, height = 8.5)
# print(fig)
# dev.off()
# Uncertainty -----
# sdr <- sdreport(obj)
# saveRDS(sdr, "data/external/mississippi_small/results/sdr.rds")
# 
# yhat_sd <- sdr$sd
# yhat_sd <- array(yhat_sd, dim = dim(data$y_n_t_v))
# 
# reaches$sd1 <- yhat_sd[, 1, 1]
# reaches$sd2 <- yhat_sd[, 2, 1]
# reaches$sd3 <- yhat_sd[, 3, 1]
# 
# 
# sdp <- lapply(c("sd1", "sd2", "sd3"), \(x) {
#   ggplot(reaches) + 
#     geom_sf(aes(color = .data[[x]])) +
#     scale_color_continuous(limits = c(.001, 1.2),
#                            breaks = c(.001, .01, .1, 1),
#                            palette = "magma", 
#                            na.value = "transparent", 
#                            transform = "log"
#     ) + 
#     theme_minimal()
# })
# fig <- ((p[[1]] + p[[2]] + p[[3]]) / (p[[4]] + p[[5]] + p[[6]]) / (sdp[[1]]+ sdp[[2]] + sdp[[3]]))
# pdf("result.pdf", width = 11, height = 11)
# print(fig)
# dev.off()
# 
