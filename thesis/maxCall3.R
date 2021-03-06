## Max call option. N = 5.

rm(list = ls())
gc()
library(rcss)

## Parameters
rate <- 0.05
d <- 0.10
vol <- 0.2
n_dec <- 10
step <- 1/3
control <- matrix(c(c(1, 1), c(2, 1)), nrow = 2, byrow = TRUE)
strike <- 100
u <- (rate - d - 0.5 * vol^2) * step

## Disturbance sampling
set.seed(1)
r_index <- matrix(c(2, 3, 4, 5, 6, 2, 3, 4, 5, 6), ncol = 2)
sigma <- vol * sqrt(step)
n_disturb <- 40000
weight <- rep(1 / n_disturb, n_disturb)
disturb <- array(0, dim = c(6, 6, n_disturb))
disturb[1,1,] <- 1
disturb[2,2,] <- rlnorm(n_disturb, u, sigma)
disturb[3,3,] <- rlnorm(n_disturb, u, sigma)
disturb[4,4,] <- rlnorm(n_disturb, u, sigma)
disturb[5,5,] <- rlnorm(n_disturb, u, sigma)
disturb[6,6,] <- rlnorm(n_disturb, u, sigma)

## Generate stochastic grid
n_path_grid <- 2000
grid_disturb <- array(0, dim = c(6, 6, n_path_grid, n_dec - 1))
grid_disturb[1,1,,] <- 1
grid_disturb[2,2,,] <- rlnorm(n_path_grid * (n_dec - 1), u, sigma)
grid_disturb[3,3,,] <- rlnorm(n_path_grid * (n_dec - 1), u, sigma)
grid_disturb[4,4,,] <- rlnorm(n_path_grid * (n_dec - 1), u, sigma)
grid_disturb[5,5,,] <- rlnorm(n_path_grid * (n_dec - 1), u, sigma)
grid_disturb[6,6,,] <- rlnorm(n_path_grid * (n_dec - 1), u, sigma)

## Generate paths
start <- c(1, 100, 100, 100, 100, 100)
path_sgrid <- PathDisturb(start, grid_disturb)
start_index <- rbind(c(1, 90, 90, 90, 90, 90),
                     c(1, 100, 100, 100, 100, 100),
                     c(1, 110, 110, 110, 110, 110))
n_grid <- 2000
grid <- StochasticGrid(path_sgrid, n_grid - nrow(start_index), 10, FALSE)
grid <- rbind(start_index, grid)

## Subgradient representation of reward
in_money <- apply(grid[,2:6], 1, max) >= strike
colIndex <- max.col(grid[,2:6]) + 1
reward <- array(0, dim = c(n_grid, 6, 2, 2, n_dec - 1))
reward[in_money,1,2,2,] <- -strike
for (gg in 1:n_grid) {
    if (in_money[gg] == TRUE) {
        reward[gg, colIndex[gg], 2,2,] <- 1
    }
}
for (tt in 1:n_dec - 1){
    reward[,,,,tt] <- exp(-rate * step * (tt - 1)) * reward[,,,,tt]
}
## Subgrad representation of scrap
scrap <- array(data = 0, dim = c(n_grid, 6, 2))
scrap[in_money,1,2] <- -strike
for (gg in 1:n_grid) {
    if (in_money[gg] == TRUE) {
        scrap[gg, colIndex[gg], 2] <- 1
    }
}
scrap <- exp(-rate * step * (n_dec - 1)) * scrap

## Bellman 
time <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, weight, r_index)
time <- proc.time() - time
value <- rowSums(bellman$value[,,2,1] * grid)

## Get primal-dual bounds
n_path <- 1500
n_subsim <- 1500
## Path disturbances 
set.seed(12345)
path_disturb <- array(0, dim = c(6, 6, n_path, n_dec - 1))
path_disturb[1, 1,,] <- 1
path_disturb[2, 2,,] <- rlnorm(n_path * (n_dec - 1), u, sigma)
path_disturb[3, 3,,] <- rlnorm(n_path * (n_dec - 1), u, sigma)
path_disturb[4, 4,,] <- rlnorm(n_path * (n_dec - 1), u, sigma)
path_disturb[5, 5,,] <- rlnorm(n_path * (n_dec - 1), u, sigma)
path_disturb[6, 6,,] <- rlnorm(n_path * (n_dec - 1), u, sigma)
## Subsim disturbances
subsim <- array(0, dim = c(6, 6, n_subsim, n_path, (n_dec - 1)))
subsim[1,1,,,] <- 1
subsim[2,2,,,] <- rlnorm(n_subsim * n_path * (n_dec - 1), u, sigma)
subsim[3,3,,,] <- rlnorm(n_subsim * n_path * (n_dec - 1), u, sigma)
subsim[4,4,,,] <- rlnorm(n_subsim * n_path * (n_dec - 1), u, sigma)
subsim[5,5,,,] <- rlnorm(n_subsim * n_path * (n_dec - 1), u, sigma)
subsim[6,6,,,] <- rlnorm(n_subsim * n_path * (n_dec - 1), u, sigma)
subsim_weight <- rep(1 / n_subsim, n_subsim)

## Reward function
RewardFunc <- function(state, time) {
    output <- array(data = 0, dim = c(nrow(state), 2, 2))
    output[,2,2] <- exp(-rate * step * (time - 1)) * pmax(apply(state[,2:6], 1, max) - 100, 0)
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(data = 0, dim = c(nrow(state), 2))
    output[,2] <- exp(-rate * step * (n_dec - 1)) * pmax(apply(state[,2:6], 1, max) - 100, 0)
    return(output)
}

## Containers
start_index <- rbind(c(1, 90, 90, 90, 90, 90),
                     c(1, 100, 100, 100, 100, 100),
                     c(1, 110, 110, 110, 110, 110))
subgradient<- rep(NA, 3)
boundsIndex <- cbind(subgradient, subgradient)

## Looping
for (ss in 1:3) {
    start <- start_index[ss,]
    path <- PathDisturb(start, path_disturb)
    policy <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
    time2 <- proc.time()
    mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
    bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, policy)
    time2 <- proc.time() - time2
    boundsIndex[ss,] <- GetBounds(bounds, 0.05, 2)
    path_nn <- rflann::FastKDNeighbour(matrix(start, nrow = 1), grid, 1)
    subgradient[ss] <- value[path_nn]
}

table <- round(cbind(subgradient, boundsIndex), 3)
print(table)
print(time)
print(time2)
