## Comparing original, neighbours and matrices as function of grid density

rm(list = ls())
gc()
library(rcss)

## Parameters
rate <- 0.06
step <- 0.02
vol <- 0.2
n_dec <- 51
control <- matrix(c(c(1, 1), c(2, 1)), nrow = 2, byrow = TRUE)
strike <- 40
grid <- as.matrix(cbind(rep(1, 301), seq(30, 60, length = 301)))

## Disturbance sampling
r_index <- matrix(c(2, 2), ncol = 2)
u <- (rate - 0.5 * vol^2) * step
sigma <- vol * sqrt(step)
condExpected <- function(a, b){
    aa <- (log(a) - (u + sigma^2)) / sigma
    bb <- (log(b) - (u + sigma^2)) / sigma
    return(exp(u + sigma^2 / 2) * (pnorm(bb) - pnorm(aa)))
}
weight <- rep(1 / 1000, 1000)
disturb <- array(0, dim = c(2, 2, 1000))
disturb[1,1,] <- 1
part <- qlnorm(seq(0, 1, length = 1000 + 1), u, sigma)
for (i in 1:1000) {
    disturb[2,2,i] <- condExpected(part[i], part[i+1]) / (plnorm(part[i+1], u, sigma) - plnorm(part[i], u, sigma))
}

## Subgradient representation of reward
in_money <- grid[,2] <= strike
reward <- array(0, dim = c(301, 2, 2, 2, n_dec - 1))       
reward[in_money,1,2,2,] <- strike
reward[in_money,2,2,2,] <- -1
for (tt in 1:n_dec - 1){
    reward[,,,,tt] <- exp(-rate * step * (tt - 1)) * reward[,,,,tt]
}
## Subgrad representation of scrap
scrap <- array(data = 0, dim = c(301, 2, 2))
scrap[in_money,1,2] <- strike
scrap[in_money,2,2] <- -1
scrap <- exp(-rate * step * (n_dec - 1)) * scrap

## Bellman 
time <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, weight, r_index)
time <- proc.time() - time
value <- rowSums(bellman$value[,,2,1] * grid)

## Get primal-dual bounds
n_path <- 1000
n_subsim <- 1000
## Path disturbances 
set.seed(12345)
path_disturb <- array(0, dim = c(2, 2, n_path, n_dec - 1))
path_disturb[1, 1,,] <- 1
rand1 <- rnorm(n_path * (n_dec - 1) / 2)
rand1 <- as.vector(rbind(rand1, -rand1))
path_disturb[2, 2,,] <- exp((rate - 0.5 * vol^2) * step + vol * sqrt(step) * rand1)
## Subsim disturbances
subsim <- array(0, dim = c(2, 2, n_subsim, n_path, (n_dec - 1)))
subsim[1,1,,,] <- 1
rand2 <- rnorm(n_subsim * n_path * (n_dec - 1) / 2)
rand2 <- as.vector(rbind(rand2, -rand2))
subsim[2,2,,,] <- exp((rate - 0.5 * vol^2) * step + vol * sqrt(step) * rand2)
subsim_weight <- rep(1 / n_subsim, n_subsim)

## Reward function
RewardFunc <- function(state, time) {
    output <- array(data = 0, dim = c(nrow(state), 2, 2))
    output[,2,2] <- exp(-rate * step * (time - 1)) * pmax(40 - state[,2], 0)
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(data = 0, dim = c(nrow(state), 2))
    output[,2] <- exp(-rate * step * (n_dec - 1)) * pmax(40 - state[,2], 0)
    return(output)
}

## Containers
start_index <- cbind(rep(1,5),c(36,38,40,42,44))
subgradient<- rep(NA,5)
primal <- cbind(subgradient, subgradient)
dual <- primal

## Looping
for (ss in 1:5) {
    start <- start_index[ss,]
    path <- PathDisturb(start, path_disturb)
    policy <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
    time2 <- proc.time()
    mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
    bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, policy)
    time2 <- proc.time() - time2
    primal[ss,1] <- mean(bounds$primal[,2,1])
    primal[ss,2] <- sd(bounds$primal[,2,1]) / sqrt(n_path)
    dual[ss,1] <- mean(bounds$dual[,2,1]) 
    dual[ss,2] <- sd(bounds$dual[,2,1]) / sqrt(n_path)
    path_nn <- rflann::FastKDNeighbour(matrix(start, nrow = 1), grid, 1)
    subgradient[ss] <- value[path_nn]
}

table <- round(cbind(subgradient, primal, dual), 5)
print(table)
print(time)
print(time2)
