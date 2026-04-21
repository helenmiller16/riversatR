library(data.table)
library(Matrix)
library(RTMB)
source(here::here("src_interpolate/functions.R"))
TIME_INTERVAL <- 8
VARIABLE <- "fai"
reaches <- sf::read_sf("data/external/mississippi_small/reaches.gpkg")
reaches <-dplyr::rename(reaches, 
                        reach_id_up = "rch_id_up", 
                        reach_id_down = "rch_id_dn")
reaches$reach_len <- reaches$reach_len / 1000 # convert to km
reflectance_dt <- fread("data/external/mississippi_small/landsat_2015_10_01-2016_10-01.csv")
reflectance_dt[, ndti := (red - green)/(red + green)]
reflectance_dt[, fai  := nir08 - (red + (swir16 - red)*(865-660)/(1600-660))]
variables = VARIABLE
graph <- make_graph(reaches)
data_ <- make_data(reflectance_dt, 
                  graph, 
                  variables = variables,
                  time_interval = TIME_INTERVAL, 
                  weight_variable = "facc")
time_step_date_map <- data_$time_step_date_map
data_$time_step_date_map <- NULL
y_n_t_v <- data_$y_n_t_v
y_n_t_v[y_n_t_v > 1 | y_n_t_v < -1] <- NA

data_$y_n_t_v <- y_n_t_v[, , ]
# data <- c(data_, 
#           log_gamma1 = 0, 
#           log_gamma2 = 0)
# params <- list(log_theta = -5, 
#                # log_gamma1 = 0, # spacetime
#                log_gamma2 = 0, # space only
#                # log_gamma3 = 0, # time only
#                log_sigma = 0, 
#                mu = 0, 
#                z_n_t = rnorm(prod(dim(data$y_n_t_v))),
#                w_n = rnorm(dim(data$y_n_t_v)[1])
#                )

# f1 <- objective_function(tailup_mean,
#                          data,
#                          params)
# obj <- MakeADFun(f1, params, random = c("w_n"))
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# r <- obj$report()
# w_n <- r$`params$w_n`

data <- c(data_)
params <- list(log_theta = -5, 
               log_gamma1 = 0, # spacetime
               log_gamma2 = 0, # space only
               # log_gamma3 = 0, # time only
               log_rhoT = -1,
               log_sigma = -3, 
               mu = 0, 
               z_n_t = rnorm(prod(dim(data$y_n_t_v))),
               w_n = rnorm(dim(data$y_n_t_v)[1])
)

f <- objective_function(tailup_iter_time,
                        data,
                        params)
# f <- objective_function(tailup_no_covariates, 
#                         data, 
#                         params)


obj <- MakeADFun(f, params, random = c("z_n_t", "w_n"), 
                 # map = list(log_sigma = factor(NA))
                 )
opt <- nlminb(obj$par, obj$fn, obj$gr)
r <- obj$report()

yhat <- r$x_n_t
yhat <- array(yhat, dim = dim(data$y_n_t_v))
plot(data$y_n_t_v, yhat)
abline(0, 1, col= "blue")

# Save yhat
preds <- as.data.table(yhat)
dates <- seq(min(reflectance_dt$date), max(reflectance_dt$date), by = TIME_INTERVAL) |> 
  as.character()

# names(preds) <- dates
preds$node_id <- 1:nrow(preds)
preds <- melt(preds, 
              id.vars = "node_id", 
              variable.name = "time_step", 
              value.name = VARIABLE)
preds[, time_step := as.numeric(sub("V", "", time_step))]
preds[, date := dates[time_step]]
nodes <- as.data.table(sfnetworks::activate(graph, "nodes"))
preds$reach_id <- nodes[.(preds$node_id), on = "node_id"]$reach_id
preds$date <- as.Date(preds$date)
preds <- merge(preds, time_step_date_map, 
               all.x = TRUE, 
               allow.cartesian=TRUE)
# for those with missing date_diff, use date_diff from closest reach
p <- unique(preds[, .(reach_id, date_diff)])[, .(date_diff = min(date_diff)), reach_id]

r_missing <- reaches[reaches$reach_id %in% p[is.na(date_diff)]$reach_id, ]
r_not_missing <- reaches[!reaches$reach_id %in% p[is.na(date_diff)]$reach_id, ]
r_not_missing$date_diff <- p[.(r_not_missing$reach_id), on = "reach_id"]$date_diff
r_missing$date_diff <- r_not_missing$date_diff[sf::st_nearest_feature(r_missing, r_not_missing)]
r <- rbind(r_missing, r_not_missing)
preds <- merge(preds, sf::st_drop_geometry(r[, c("reach_id", "date_diff")]))
preds[, date_diff := date_diff.x]
preds[is.na(date_diff), date_diff := date_diff.y]
preds$date_diff.x <- NULL
preds$date_diff.y <- NULL

fwrite(preds, "data/external/mississippi_small/results/yhat_iter_time_fai.csv")

# bad_alloc here

# https://groups.google.com/g/tmb-users/c/9mlEeG_D430/m/EROcIHTaCAAJ
h <- obj$env$spHess()
sdr <- sdreport(obj)
saveRDS(summary(sdr, "random"), paste0("data/external/mississippi_small/results/sdr_random_iter_time.rds"))
saveRDS(summary(sdr, "fixed"), paste0("data/external/mississippi_small/results/sdr_fixed_iter_time.rds"))

yhat1 <- readRDS("data/external/mississippi_small/results/yhat_1.rds")
yhat2 <- readRDS("data/external/mississippi_small/results/yhat_4.rds")
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
# obj <- MakeADFun(f, params, random = "z_n_t"
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# saveRDS(opt, "data/external/mississippi_small/results/opt.rds")
# r <- obj$report()
# yhat <- r$x_n_t
# yhat <- array(yhat, dim = dim(data$y_n_t_v))
# saveRDS(yhat, "data/external/mississippi_small/results/yhat.rds")
# plot(data$y_n_t_v[, , 1], yhat)
# abline(0,1, col = "blue")
reaches$obs1 <- data$y_n_t_v[, 7]
reaches$obs2 <- data$y_n_t_v[, 8]
reaches$obs3 <- data$y_n_t_v[, 9]
reaches$yhat1 <- yhat[, 7]
reaches$yhat2 <- yhat[, 8]
reaches$yhat3 <- yhat[, 9]
reaches$mean  <- r$w_n

library(ggplot2)
library(patchwork)
p <- lapply(c("obs1", "obs2", "obs3",
              "yhat1", "yhat2", "yhat3"), \(x) {
                ggplot(reaches) +
                  geom_sf(aes(color = .data[[x]])) +
                  scale_color_continuous(limits = c(-.5, .5),
                                         palette = "viridis",
                                         na.value = "transparent") +
                  theme_minimal()
              })
# 
# 
fig <- (p[[1]] + p[[2]] + p[[3]]) / (p[[4]] + p[[5]] + p[[6]])
pdf("fitted_7_itertime4.pdf", width = 11, height = 8.5)
print(fig)
dev.off()

# Look at a few time series

plot(yhat[12,], type = "l")
points(data$y_n_t_v[11,])


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
