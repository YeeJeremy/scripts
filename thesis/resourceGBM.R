## Mine valuation under GBM

rm(list = ls())
gc()
library(rcss)

## CSS Algorithm
n_grid <- 2001
grid <- as.matrix(cbind(rep(1, n_grid), seq(0, 10, length = n_grid)))
rate <- 0.1
delta <- 0.01
vol <- sqrt(0.08)
step <- 0.25

## Disturbance sampling
n_disturb <- 1000
disturb_weight <- rep(1/n_disturb, n_disturb)
disturb <- array(0, dim = c(2, 2, n_disturb))
disturb[1, 1,] <- 1
quantile <- rep(NA, n_disturb)
u <- (rate - delta - 0.5 * vol^2) * step
sigma <- vol * sqrt(step)
part <- qlnorm(seq(0, 1, length = n_disturb + 1), u, sigma)
condExpected <- function(a, b){  ##[a,b]
    aa <- (log(a) - (u+sigma^2))/sigma
    bb <- (log(b) - (u+sigma^2))/sigma
    vv <- exp(u + sigma^2/2) * (pnorm(bb) - pnorm(aa))
    return(vv)
}
for (i in 1:n_disturb) {
    quantile[i] <- condExpected(part[i], part[i+1]) / (plnorm(part[i+1], u, sigma) - plnorm(part[i], u, sigma))
}
disturb[2, 2,] <- quantile
r_index <- matrix(c(2, 2), ncol = 2)

## Control matrix, a = 1 (close), 2 (abandon), 3 (open)
## p = 1 (exhausted), levels + 1 (full open), 2 * levels + 1 (full closed)
levels <- 15 / step
control <- matrix(data = 1, nrow = (2 * levels + 1), ncol = 3)
control[2:(2 * levels + 1), 1] <- (levels + 2):(2 * levels + 1)
control[2:(2 * levels + 1), 3] <- 1:levels

## Subgrad rep of rewards
H1 <- 5 * step
H2 <- -2.5 * step
maint <- 0.5 * step
switch <- 0.2
inflation <- 0.08
tax_property <- 0.02
n_dec <- 1 + 30 / step + 1
discount <- exp(-(rate + tax_property) * step)
n_pos <- nrow(control)
n_action <- ncol(control)
n_dim <- ncol(grid)
reward <- array(data = 0, dim = c(n_grid, n_dim, n_action, n_pos, n_dec - 1))
for (p in 2:n_pos) {
    for (t in 1:(n_dec - 1)) {
        pi <- exp(inflation * (t-1) * step)
        adjust <- discount ^ (t-1)
        ## Close the asset
        if (p > (levels + 1)) {  ## Closed
            ## Close
            reward[, 1, 1, p, t] <- adjust * -maint * pi
            ## Open
            reward[, 1, 3, p, t] <- (H2 - switch) * pi * adjust
            reward[, 2, 3, p, t] <- H1 * adjust
        } else if (p <= (levels + 1)) {  ## Opened
            ## Close
            reward[, 1, 1, p, t] <- -(maint + switch) * pi * adjust
            ## Open
            reward[, 1, 3, p, t] <- H2 * pi * adjust
            reward[, 2, 3, p, t] <- H1 * adjust
        }
    }
}

## Subgrad rep of scrap
scrap <- array(data = 0, dim = c(n_grid, n_dim, n_pos))

## Performing fast bellman recursion
time1 <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, disturb_weight, r_index)
time1 <- proc.time() - time1

## Reward function
RewardFunc <- function(state, time) {
    Fpi <- exp(inflation * (time-1) * step)
    Fdiscount <- exp(-(rate + tax_property) * (time-1) * step)
    FH1 <- 5 * step * Fdiscount
    FH2 <- -2.5 * step * Fdiscount * Fpi
    Fmaint <- 0.5 * step *  Fdiscount * Fpi
    Fswitch <- 0.2 * Fdiscount * Fpi
    output <- array(0, dim = c(nrow(state), n_action, n_pos))
    output[,1,(2:(levels+1))] <- -(Fmaint + Fswitch)
    output[,3,(2:(levels+1))] <- FH2 + FH1 * state[,2] 
    output[,1,(levels+2):n_pos] <- -Fmaint
    output[,3,(levels+2):n_pos] <- FH2 - Fswitch + FH1 * state[,2]       
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(data = 0, dim = c(nrow(state), n_pos))
    return(output)
}

### Perform the solution diagnostics
## Path disturbances
set.seed(12345)
n_path <- 500
path_disturb <- array(0, dim = c(2, 2, n_path, n_dec - 1))
path_disturb[1,1,,] <- 1
rand1 <- rnorm(n_path * (n_dec - 1) / 2)
rand1 <- as.vector(rbind(rand1, -rand1))
path_disturb[2,2,,] <- exp((rate - delta - 0.5 * vol^2) * step + vol * sqrt(step) * rand1)

## Subsimulation disturbance
n_subsim <- 500
subsim_weight <- rep(1 / n_subsim, n_subsim)
subsim <- array(0, dim = c(2, 2, n_subsim, n_path, n_dec - 1))
subsim[1,1,,,] <- 1
rand2 <- rnorm(n_subsim * n_path * (n_dec - 1) / 2)
rand2 <- as.vector(rbind(rand2, -rand2))
subsim[2,2,,,] <- exp((rate - delta - 0.5 * vol^2) * step + vol * sqrt(step) * rand2)

## Looping
startIndex <- cbind(rep(1,8), seq(0.3, 1, by = 0.1))
pp <- levels + 1
results <- matrix(NA, nrow=nrow(startIndex), ncol = 6)
for (i in 1:nrow(startIndex)) {
    ## Paths
    start <- startIndex[i,]
    path <- PathDisturb(start, path_disturb)
    policy <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
    ## Primal and dual values
    time2 <- proc.time()
    mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
    bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, policy)
    time2 <- proc.time() - time2
    ## Store
    ind <- rflann::FastKDNeighbour(matrix(start, nrow = 1), grid, 1)
    results[i,1] <- sum(bellman$value[ind, , pp, 1] * grid[ind,])
    results[i,4] <- sum(bellman$value[ind, , n_pos, 1] * grid[ind,])
    results[i,2:3] <- GetBounds(bounds, 0.01, pp)
    results[i,5:6] <- GetBounds(bounds, 0.01, n_pos)
}

## Print results
print(round(results,3))
print(time1)
print(time2)
