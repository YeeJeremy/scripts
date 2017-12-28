## Mine valuation under AR process

rm(list = ls())
gc()
library(rcss)

## Parameters
n_grid <- 2000
grid <- as.matrix(cbind(rep(1, n_grid), seq(-5, 5, length = n_grid)))
rate <- 0.1
delta <- 0.01
step <- 0.25
vol <- sqrt(0.08)
n_dec <- 1 + 30 / step + 1
mu <- (rate - delta - 0.5 * vol^2) * step 
psi <- 0.6

## Disturbances
n_disturb <- 1000
disturb_weight <- rep(1/n_disturb, n_disturb)
disturb <- array(matrix(c(1, 0, 0, psi), ncol = 2, byrow = TRUE), dim = c(2, 2, n_disturb))
part <- qnorm(seq(0, 1, length = n_disturb + 1))
condExpected <- function(a, b){ ## (a,b)
    return(1/sqrt(2 * pi) * (exp(-a^2/2)- exp(-b^2/2)))
}
for (i in 1:n_disturb) {
    disturb[2,1,i] <- mu + vol * sqrt(step) *condExpected(part[i], part[i+1]) / (pnorm(part[i+1]) - pnorm(part[i]))
}
r_index <- matrix(c(2, 1), ncol = 2)

## Control matrix
levels <- 15 / step
control <- array(data = 1, dim = c(2 * levels + 1, 3))
control[2:(2 * levels + 1), 1] <- (levels + 2):(2 * levels + 1)
control[2:(2 * levels + 1), 3] <- 1:levels

## Rewards array
H1 <- 5 * step
H2 <- -2.5 * step
maint <- 0.5 * step
switch <- 0.2
inflation <- 0.08
tax_property <- 0.02
discount <- exp(-(rate + tax_property) * step)
n_pos <- nrow(control)
n_action <- ncol(control)
n_dim <- ncol(grid)
reward <- array(0, dim = c(n_grid, n_dim, n_action, n_pos, n_dec - 1))
slope <- H1 * exp(grid[,2])
## Assigning values to reward array
for (p in 2:n_pos) {
    for (t in 1:(n_dec-1)) {
        pi <- exp(inflation * (t-1) * step)
        adjust <- discount^(t-1)
        intercept1 <- (slope + H2 * pi) * adjust - (slope * adjust * grid[,2])
        intercept2 <- (slope + (H2 - switch) * pi) * adjust - (slope * adjust * grid[,2])
        if (p > (levels + 1)) {  ## Closed
            ## Close
            reward[, 1, 1, p, t] <- adjust * -maint * pi
            ## Open
            reward[, 1, 3, p, t] <- intercept2
            reward[, 2, 3, p, t] <- slope * adjust
        } else if (p <= (levels + 1)) {  ## Opened
            ## Close
            reward[, 1, 1, p, t] <- -(maint + switch) * pi * adjust
            ## Open
            reward[, 1, 3, p, t] <- intercept1
            reward[, 2, 3, p, t] <- slope * adjust
        }
    }
}

## Subgrad rep of scrap
scrap <- array(data = 0, dim = c(n_grid, n_dim, n_pos))

## Performing fast bellman recursion
time1 <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, disturb_weight, r_index)
time1 <- proc.time() - time1

## Solution diagnostics
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
    output[,3,(2:(levels+1))] <- FH2 + FH1 * exp(state[,2]) 
    output[,1,(levels+2):n_pos] <- -Fmaint
    output[,3,(levels+2):n_pos] <- FH2 - Fswitch + FH1 * exp(state[,2])       
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(data = 0, dim = c(nrow(state), n_pos))
    return(output)
}

## Path and subsim disturbances
n_path <- 500
path_disturb <- array(matrix(c(1, 0, 0, psi), ncol = 2, byrow = TRUE), dim = c(2, 2, n_path, n_dec - 1))
n_subsim <- 500
subsim_weight <- rep(1 / n_subsim, n_subsim)
subsim <- array(matrix(c(1, 0, 0, psi), ncol = 2, byrow = TRUE), dim = c(2, 2, n_subsim, n_path, n_dec - 1))

## Looping
startIndex <- c(0.3, 0.4, 0.5, 0.6, 0.7)
pp <- levels + 1
results <- matrix(NA, nrow = length(startIndex), ncol = 6)
ind <- rflann::FastKDNeighbour(matrix(log(startIndex), ncol = 1), as.matrix(grid[,2]), 1)

for (z in 1:length(startIndex)) {
    ## Generate and classify the paths
    set.seed(12345)
    rand1 <- rnorm(n_path * (n_dec - 1) / 2)
    rand1 <- as.vector(rbind(rand1, -rand1))  ## anti-thetic disturbances
    path_disturb[2,2,,] <- mu + rand1 * vol * sqrt(step)
    start <- c(1, log(startIndex[z]))
    path <- PathDisturb(start, path_disturb)
    policy <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
    ## Specifying subsimulation disturbance matrices    
    rand2 <- rnorm(n_subsim * n_path * (n_dec - 1)/ 2)
    rand2 <- as.vector(rbind(rand2, -rand2))
    subsim[2,1,,,] <- mu + rand2 * vol * sqrt(step)
    ## Primal and dual values
    time2 <- proc.time()
    mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
    bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, policy)
    time2 <- proc.time() - time2
    ## Storing the results
    results[z,1] <- sum(bellman$value[ind[z], , pp, 1] * grid[ind[z],])
    results[z,4] <- sum(bellman$value[ind[z], , n_pos, 1] * grid[ind[z],])
    results[z,2:3] <- GetBounds(bounds, 0.01, pp)
    results[z,5:6] <- GetBounds(bounds, 0.01, n_pos)
}

## Print results
print(round(results,3))
print(time1)
print(time2)
