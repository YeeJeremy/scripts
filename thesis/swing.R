## Swing option

rm(list = ls())
gc()
library(rcss)

## Set Parameters
rho <- 0
kappa <- 0.9
mu <- 0
sigma <- 0.5
K <- 0 
n_dec <- 1001  ## number of time epochs
N <- 100  ## number of rights
n_pos <- N + 1 ## number of positions

## Grid
n_grid <- 101
grid <- cbind(rep(1, n_grid), seq(-2, 2, length = n_grid))

## Control matrix
control <- cbind(c(1, 1:N), 1:(N + 1))

## Reward subgradient representation
reward <- array(0, dim = c(n_grid, 2, 2, nrow(control), n_dec - 1))
slope <- exp(grid[, 2])
for (tt in 1:(n_dec - 1)) {
    discount <- exp(-rho * (tt - 1))
    for (pp in 2:n_pos) {
        intercept <- (exp(grid[,2]) - K * discount) - slope * grid[, 2]
        reward[, 1, 1, pp, tt] <- intercept
        reward[, 2, 1, pp, tt] <- slope
    }
}
## Scrap subgradient representation
scrap <- array(0, dim = c(n_grid, 2, nrow(control)))
discount <- exp(-rho * (n_dec - 1))
for (pp in 2:n_pos) {
    intercept <- (exp(grid[,2]) - K * discount) - slope * grid[, 2]
    scrap[, 1, pp] <- intercept
    scrap[, 2, pp] <- slope
}

## Disturbance sampling
n_disturb <- 1000
weight <- rep(1/n_disturb, n_disturb)
disturb <- array(0, dim = c(2, 2, n_disturb))
disturb[1, 1,] <- 1
disturb[2, 2,] <- 1 - kappa
CondExpected <- function(a, b){
    return(1/sqrt(2 * pi) * (exp(-a^2/2)- exp(-b^2/2)))
}
part <- qnorm(seq(0, 1, length = n_disturb + 1))
for (i in 1:n_disturb) {
    disturb[2,1,i] <- kappa * mu + sigma * (CondExpected(part[i], part[i+1]) / (pnorm(part[i+1]) - pnorm(part[i])))
}

## Bellman recursion
r_index <- matrix(c(2, 1), ncol = 2)
time1 <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, weight, r_index)
time1 <- proc.time() - time1

## Exact reward function
RewardFunc <- function(state, time) {
    output <- array(0, dim = c(nrow(state), 2, nrow(control)))
    discount <- exp(-rho * (time - 1))
    for (i in 2:nrow(control)) {
        output[, 1, i] <- pmax(exp(state[, 2]) - K * discount, 0)
    }
    return(output)
}
## Exact scrap function
ScrapFunc <- function(state) {
    output <- array(0, dim = c(nrow(state), nrow(control)))
    discount <- exp(-rho * (n_dec - 1))
    for (i in 2:nrow(control)) {
        output[, i] <- pmax(exp(state[, 2]) - K * discount, 0)
    }
    return(output)
}

## Solution diagnostics
set.seed(12345)
## Generate paths
n_path <- 100
path_disturb <- array(0, dim = c(2, 2, n_path, n_dec - 1))
path_disturb[1, 1,,] <- 1
path_disturb[2, 2,,] <- 1 - kappa
rand1 <- rnorm(n_path * (n_dec - 1) / 2)
rand1 <- as.vector(rbind(rand1, -rand1))
path_disturb[2, 1,,] <- kappa * mu + sigma * rand1
start <- c(1, 0)
path <- PathDisturb(start, path_disturb)
policy <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
## Set subsimulation disturbances
n_subsim <- 100
subsim <- array(0, dim = c(2, 2, n_subsim, n_path, n_dec - 1))
subsim[1, 1,,,] <- 1
subsim[2, 2,,,] <- 1 - kappa
rand2 <- rnorm(n_subsim * n_path * (n_dec - 1) / 2)
rand2 <- as.vector(rbind(rand2, -rand2))
subsim[2, 1,,,] <- kappa * mu + sigma * rand2
subsim_weight <- rep(1 / n_subsim, n_subsim)

## Primal-dual
time2 <- proc.time()
mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, policy)
time2 <- proc.time() - time2


index <- c(2:6,11,16,seq(21,101, by= 10))
output <- matrix(NA, ncol = 2, nrow = length(index))
value <- rep(NA, length(index))
valueF <- matrix(NA, ncol = length(index), nrow = n_grid)
for (i in 1:length(index)) {
    output[i,] <- GetBounds(bounds, 0.01, index[i])
    value[i] <- sum(bellman$value[51,,index[i],1] * grid[51,])
    valueF[,i] <- rowSums(bellman$value[,,index[i],1] * grid)
}

print(round(cbind(index-1,value,output), 3))
print(time1)
print(time2)

##setEPS()
##postscript("Swing1.eps")
matplot(exp(grid[,2]),valueF[,1:5],type="l", xlab="S_0", ylab ="Value")
##dev.off()

##setEPS()
##postscript("Swing2.eps")
matplot(exp(grid[,2]),valueF[,11:16],type="l", xlab="S_0", ylab ="Value")
##dev.off()
