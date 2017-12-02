# test data (True gaussian mixture)
library(mvtnorm)
# library(pldensity)
# set.seed(110104)
nclust <- 4
d <- 2
mu <- list(c(.75, .5), c(.5, .75), c(.25, .5), c(.75, .3))
S <- list(matrix(.1^2 * c(1, .5, .5, 1), 2, 2), .05^2 * diag(2),
          matrix(.1^2 * c(1, -.5, -.5, 1), 2, 2), .075^2 * diag(2))
n <- 1000
k <- sample(1:4, n, TRUE)
y <- matrix(0, n, d)
for (i in 1:nclust) {
  idx <- (k == i)
  print(mu[[i]])
  y[idx, ] <- rmvnorm(sum(idx), mu[[i]], S[[i]])
}
summary(y)
plot(y, xlim = c(0, 1), ylim = c(0, 1), bg = k, pch = 21,
     main = "Simulated normal mixture data")

system.time({
  res2 <- dp_normal_mix(
    y[ , ], 
    N = 500,
    alpha = 10, 
    lambda = runif(2), 
    kappa = 1, 
    nu = 2,
    Omega =  0.1 ^ 2 * diag(2))  
})


res <- res2[[1]]
points(t(res$mu), bg = "yellow", pch = 21)

resol <- 100
mesh <- expand.grid(x = seq(0, 1, length.out = resol), y = seq(0, 1, length.out = resol))
mesh$z <- 0

m <- sum(res$n[res$n > 1])
idx <- which(res$n > 1)
for (i in idx) {
  # if (res$n[i] > 1) {
  #   covmat <-  res$S[ , ,i] / (res$n[i] - 2) 
  # } else {
  #   covmat <- 1 * 0.15 ^ 2 * diag(2)
  # }
  mesh$z <- mesh$z + (res$n[i] / res$m) * dmvnorm(mesh[ ,c("x", "y")], res$mu[ , i, drop = TRUE], res$S[ , ,i] / (res$n[i] - 2))
}
mesh$z <- mesh$z / max(mesh$z)
z <- matrix(mesh$z, resol, resol)
contour(z, add = TRUE)
# 
library(plotly)
plot_ly(z = ~t(z)) %>% add_heatmap()
# 
# plot_ly(z = ~t(z)) %>% add_surface()


# =================================

# 0. Libraries ===================================
library(tidyverse)
library(lubridate)

# 1. Read Data =====================================
rides <- read_csv("C:/Users/mbg877/Google Drive/P1_Ride_Austin/00_Data/Clean_Database/Rides_A.csv") 
rides <- rides %>%
  mutate(datehour = ymd_h(paste(paste(year(started_on), month(started_on), day(started_on), sep = "-"), hour(started_on))))

rides_count <- rides %>% 
  group_by(datehour) %>% 
  count()

plot(tail(rides_count$datehour, 100), tail(rides_count$n[6800:7027], 100), type = "l")
abline(h = mean(rides_count$n), col = "red")


example <- rides %>% 
  filter(datehour == ymd_h("2017-04-06 24"))

x <- example %>% 
  select(start_location_lat, start_location_long) %>% 
  rename(lat = start_location_lat, lon = start_location_long)

x$lat <- 0.1 + 0.8 * (x$lat - min(x$lat)) / (max(x$lat) - min(x$lat))
x$lon <- 0.1 + 0.8 * (x$lon - min(x$lon)) / (max(x$lon) - min(x$lon))
x <- data.matrix(x)

plot(x, xlim = c(0, 1), ylim = c(0, 1), main = "Ride Austin Trips")

system.time({
  res2 <- dp_normal_mix(
    x, 
    N = 500,
    alpha = 50, 
    lambda = runif(2), 
    kappa = 1, 
    nu = 2,
    Omega =  0.1 ^ 2 * diag(2))  
})

res <- res2[[1]]
points(t(res$mu), bg = "yellow", pch = 21)
m <- sum(res$n[res$n > 1])
idx <- which(res$n > 1)
for (i in idx) {
  # if (res$n[i] > 1) {
  #   covmat <-  res$S[ , ,i] / (res$n[i] - 2) 
  # } else {
  #   covmat <- 1 * 0.15 ^ 2 * diag(2)
  # }
  mesh$z <- mesh$z + (res$n[i] / res$m) * dmvnorm(mesh[ ,c("x", "y")], res$mu[ , i, drop = TRUE], res$S[ , ,i] / (res$n[i] - 2))
}
mesh$z <- mesh$z / max(mesh$z)
z <- matrix(mesh$z, resol, resol)
contour(z, add = TRUE)
# 
library(plotly)
plot_ly(z = ~t(z)) %>% add_heatmap()

