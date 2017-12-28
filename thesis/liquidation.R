## Optimal liquidation


rm(list = ls())
gc()
library(rcss)
set.seed(12345)

## Parameters for bid-ask spread process
phi <- 0.9
sigma <- 0.5
mu <- 1
n_dec <- 51
n_grid <- 200 + 1
grid <- cbind(rep(1, n_grid), seq(5, -5, length = n_grid))

## Disturbance samplings
n_disturb <- 1000
weight <- rep(1/n_disturb, n_disturb)
disturb <- array(0, dim = c(2, 2, n_disturb))
disturb[1,1,] <- 1
disturb[2,2,] <- phi
CondExpected <- function(a, b){
    return(1/sqrt(2 * pi) * (exp(-a^2/2)- exp(-b^2/2)))
}
part <- qnorm(seq(0, 1, length = n_disturb + 1))
for (i in 1:n_disturb) {
    disturb[2,1,i] <- sigma * (CondExpected(part[i], part[i+1]) / (pnorm(part[i+1]) - pnorm(part[i])))
}
r_index <- matrix(c(2, 1), ncol = 2)


## Rewards subgradient arrays
n_pos <- 51  ## number of assets in portfolio. pos = 1 --> no asset left
max_order <- 3  ## max order that can be placed
n_action <- 1 + 2 * max_order  ## number of actions, action = 1 is do nothing
benchmark <- seq(n_pos, 1, length = n_dec)  ## benchmark across decision epoch
benchmark_range <- cbind(pmax(1, benchmark - 1), benchmark + 1)  ## benchmark range
benchmark_range[n_dec,] <- 1 ## Must not be any assets left
penalty_loop <- rbind(seq(1, 1, length = n_dec), seq(0, 1, length = n_dec), seq(0, 2.5, length = n_dec))
reward1 <- array(0, dim = c(n_grid, 2, n_action, n_pos, n_dec - 1))
reward2 <- array(0, dim = c(n_grid, 2, n_action, n_pos, n_dec - 1))
scrap1 <- array(0, dim = c(n_grid, 2, n_pos))
scrap2 <- array(0, dim = c(n_grid, 2, n_pos))

## Control matrix
## [i,j,k]: i = current pos, j = action, k = prob of being in that pos next
control <- array(0, dim = c(n_pos, n_action, n_pos))
control[1,,1] <- 1  ## absorbing state (no more assets to liquidate)
control[,1,] <- diag(n_pos) ## Probabilities associated with doing nothing
## Assigning prob for market order
success_market <- c(1, 0.9, 0.8)
for (i in 1:max_order) {
    for (j in 2:n_pos) {
        control[j,i + 1,max(1, j - i)] <- success_market[i]  ## transaction success
        control[j,i + 1,j] <- 1 - success_market[i]  ## no transaction
    }
}
## Assigning prob for limit orders
success_limit <- c(0.3, 0.2, 0.1)  ## probability of transaction for limit order
for (i in (max_order + 1):(2 * max_order)){
    for (j in 2:n_pos) {
        ## successful transaction
        control[j,i + 1,max(1, j - (i - max_order))] <- success_limit[i - max_order]
        control[j,i + 1,j] <- 1 - success_limit[i - max_order]  ## unsuccessful
    }
}

## Sample disturbances and subsimulations for diagnostic checking
## Generate the sample disturbances
n_path <- 200
path_disturb <- array(data = matrix(c(1, 0, 0, phi), ncol = 2, byrow = TRUE),
                      dim = c(2, 2, n_path, (n_dec - 1)))
rand1 <- rnorm(n_path * (n_dec - 1) / 2)
rand1 <- as.vector(rbind(rand1, -rand1))
path_disturb[2,1,,] <- sigma * rand1
## Generate the subsimulation disturbances
n_subsim <- 200
subsim_disturb <- array(data = matrix(c(1, 0, 0, phi), ncol = 2, byrow = TRUE),
                        dim = c(2, 2, n_subsim, n_path, n_dec - 1))
rand2 <- rnorm(n_subsim * n_path * (n_dec - 1) / 2)
rand2 <- as.vector(rbind(rand2, -rand2))
subsim_disturb[2,1,,,] <- sigma * rand2
subsim_weight <- rep(1/n_subsim, n_subsim)
start_loop <- matrix(c(1, 1, 1, -1, 0, 1), ncol = 2)

## Storage containers
value1 <- matrix(NA, nrow = n_grid, ncol = nrow(penalty_loop))
value2 <- matrix(NA, nrow = n_grid, ncol = nrow(penalty_loop))
lower1 <- matrix(NA, nrow = nrow(start_loop), ncol = nrow(penalty_loop))
upper1 <- lower1
lower2 <- lower1
upper2 <- lower1
alpha <- 0.01

## Looping through the various parameter values
for (i in 1:nrow(penalty_loop)) {
    penalty <- penalty_loop[i,]
    ## Reward for single target benchmark
    ## Scrap reward1
    for (p in 2:n_pos) {
        ## Penalty from deviating from benchmark
        scrap1[,1,p] <- -penalty[n_dec] * abs(benchmark[n_dec] - p) 
    }
    ## Rewards for other times
    for (t in 1:(n_dec - 1)) {
        for (p in 2:n_pos) {
            ## Penalty from deviating from benchmark
            reward1[,1,,p,t] <- -penalty[t] * abs(benchmark[t] - p) 
        }
        ## Penalty from market orders
        reward1[,2,2:(max_order+1),2:n_pos,t] <- -1
        reward1[,1,2:(max_order+1),2:n_pos,t] <- reward1[,1,2:(max_order+1),2:n_pos,t] - mu
    }
    ## Reward for ranged benchmark
    ## Scrap reward2
    for (p in 2:n_pos) {
        ## Penalty from deviating from benchmark range
        if ((p < benchmark_range[n_dec,1]) || (p > benchmark_range[n_dec,2])) {
            scrap2[,1,p] <- -penalty[n_dec] * abs(benchmark[n_dec] - p)
        }
    }
    ## Rewards for other times
    for (t in 1:(n_dec - 1)) {
        for (p in 2:n_pos) {
            ## Penalty from deviating from benchmark range
            if ((p < benchmark_range[t,1]) || (p > benchmark_range[t,2])) {
                reward2[,1,,p,t] <- -penalty[t] * abs(benchmark[t] - p)
            }
        }
        ## Penalty from market orders
        reward2[,2,2:(max_order+1),2:n_pos,t] <- -1
        reward2[,1,2:(max_order+1),2:n_pos,t] <- reward2[,1,(max_order + 2):n_action,2:n_pos,t] - mu
    }
    ## Bellman recursions
    bellman1 <- FastBellman(grid, reward1, scrap1, control, disturb, weight, r_index)
    time1 <- proc.time()
    bellman2 <- FastBellman(grid, reward2, scrap2, control, disturb, weight, r_index)
    time1 <- proc.time() - time1
    ## Storing the values
    value1[,i] <- rowSums(bellman1$value[,,n_pos,1] * grid)
    value2[,i] <- rowSums(bellman2$value[,,n_pos,1] * grid)
    ## Exact reward functions
    ## The reward function for single target
    RewardFunc1 <- function(state, time) {
        output <- array(data = 0, dim = c(nrow(state), n_action, n_pos))
        for (p in 2:n_pos) {
            output[,,p] <- -penalty[time] * abs(benchmark[time] - p)
            output[,2:(max_order + 1), p] <- output[,2:(max_order + 1), p] - state[,2] - mu  
        }
        return(output)
    }
    ScrapFunc1 <- function(state) {
        output <- array(data = 0, dim = c(nrow(state), n_pos))
        for (p in 2:n_pos) {
            output[,p] <- -penalty[n_dec] * abs(benchmark[n_dec] - p) 
        }
        return(output)
    }
    ## The reward function for ranged target
    RewardFunc2 <- function(state, time) {
        output <- array(data = 0, dim = c(nrow(state), n_action, n_pos))
        for (p in 2:n_pos) {
            if ((p < benchmark_range[time,1]) || (p > benchmark_range[time,2])) {
                output[,,p] <- -penalty[time] * abs(benchmark[time] - p)
            }
            output[,2:(max_order + 1),p] <- output[,2:(max_order + 1),p] - state[,2] - mu
        }
        return(output)
    }
    ScrapFunc2 <- function(state) {
        output <- array(data = 0, dim = c(nrow(state), n_pos))
        for (p in 2:n_pos) {
            if ((p < benchmark_range[n_dec,1]) || (p > benchmark_range[n_dec,2])) {
                output[,p] <- -penalty[n_dec] * abs(benchmark[n_dec] - p)
            }
        }
        return(output)
    }
    ## Diagnostic checking
    for (j in 1:nrow(start_loop)) {
        ## Generate the sample paths
        start <- start_loop[j,]
        path <- PathDisturb(start, path_disturb)
        policy1 <- FastPathPolicy(path, grid, control, RewardFunc1, bellman1$expected)
        mart1 <- FastAddDual(path, subsim_disturb, subsim_weight, grid, bellman1$value, ScrapFunc1)
        duality1 <- AddDualBounds(path, control, RewardFunc1, ScrapFunc1, mart1, policy1)
        time2 <- proc.time()
        policy2 <- FastPathPolicy(path, grid, control, RewardFunc2, bellman2$expected)
        mart2 <- FastAddDual(path, subsim_disturb, subsim_weight, grid, bellman2$value, ScrapFunc2)
        duality2 <- AddDualBounds(path, control, RewardFunc2, ScrapFunc2, mart2, policy2)
        time2 <- proc.time() - time2
        ## Storing the duality results
        lower1[j,i] <- GetBounds(duality1, alpha, n_pos)[1]
        upper1[j,i] <- GetBounds(duality1, alpha, n_pos)[2]
        lower2[j,i] <- GetBounds(duality2, alpha, n_pos)[1]
        upper2[j,i] <- GetBounds(duality2, alpha, n_pos)[2]
    }
}

## Tabling results
results <- matrix(NA, nrow = (nrow(start_loop) * nrow(penalty_loop)), ncol = 6)
## CSS values
counter <- 1
for (i in 1:nrow(penalty_loop)) {
    for (j in 1:nrow(start_loop)) {
        results[counter, 1] <- value1[grid[,2] == start_loop[j,2], i]
        results[counter, 4] <- value2[grid[,2] == start_loop[j,2], i]
        counter <- counter + 1
    }
}
results[,2:3] <- cbind(c(lower1), c(upper1))
results[,5:6] <- cbind(c(lower2), c(upper2))
print(round(results, 3))
print(time1)
print(time2)

