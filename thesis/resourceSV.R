## Mine valuation under stochastic volatility

rm(list = ls())
gc()
library(rcss)
set.seed(12345)

## Parameters
rate <- 0.1
step <- 0.25
n_dec <- 30 / step + 2
beta1 <- 0.8
beta2 <- 0.1
alpha <- sqrt(0.08) * (1 - beta1 - beta2)
mu <- 0.05
start <- c(1, sqrt(0.08), 1, log(0.5))
phi <- 0.8
start_index <- rbind(c(1, sqrt(0.08), 1, log(0.3)),
                     c(1, sqrt(0.08), 1, log(0.4)),
                     c(1, sqrt(0.08), 1, log(0.5)),
                     c(1, sqrt(0.08), 1, log(0.6)),
                     c(1, sqrt(0.08), 1, log(0.7)))

## Stochastic grid parameters
n_dim <- length(start)
n_path_grid <- 2000
grid_disturb <- array(data = matrix(c(1,         0,         0,   0,
                                      alpha,     beta1, beta2,   0,
                                      0,         0,         0,   0,
                                      0,         0,         0,   phi),
                                    ncol = 4, byrow = TRUE),
                      dim = c(n_dim, n_dim, n_path_grid, (n_dec - 1)))
noise <- rnorm(n_path_grid * (n_dec - 1))
for (i in 1:(n_dec - 1)) {
    grid_disturb[3, 1,,] <- alpha * noise^2
    grid_disturb[3, 2,,] <- beta1 * noise^2
    grid_disturb[3, 3,,] <- beta2 * noise^2    
    grid_disturb[4, 1,,] <- mu * step + alpha * noise * sqrt(step)
    grid_disturb[4, 2,,] <- beta1 * noise * sqrt(step)
    grid_disturb[4, 3,,] <- beta2 * noise * sqrt(step)
}

## Generate stochastic grid
n_grid <- 3000
path_sgrid <- PathDisturb(start, grid_disturb)
grid <- StochasticGrid(path_sgrid, n_grid - 4, 10, TRUE)
grid <- rbind(start_index[1,],
              start_index[2,],
              start_index[4,],
              start_index[5,],
              grid)

## Disturbances
n_disturb <- 10000
disturb_weight <- rep(1/n_disturb, n_disturb)
disturb <- array(data = grid_disturb[,,1,1], dim = c(4, 4, n_disturb))
quantile <- rep(NA, n_disturb)
part <- qnorm(seq(0, 1, length = n_disturb + 1))
condExpected <- function(a, b){ ## (a,b)
    return(1/sqrt(2 * pi) * (exp(-a^2/2)- exp(-b^2/2)))
}
for (i in 1:n_disturb) {
    quantile[i] <- condExpected(part[i], part[i+1]) / (pnorm(part[i+1]) - pnorm(part[i]))
}
disturb[3, 1,] <- alpha * quantile^2
disturb[3, 2,] <- beta1 * quantile^2
disturb[3, 3,] <- beta2 * quantile^2
disturb[4, 1,] <- mu * step + alpha * quantile * sqrt(step)
disturb[4, 2,] <- beta1 * quantile * sqrt(step)
disturb[4, 3,] <- beta2 * quantile * sqrt(step)
r_index <- matrix(c(3, 1,
                    3, 2,
                    3, 3,
                    4, 1,
                    4, 2,
                    4, 3),
                  ncol = 2, byrow = TRUE)

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
n_position <- nrow(control)
n_action <- ncol(control)
reward <- array(0, dim = c(n_grid, n_dim, n_action, n_position, n_dec - 1))
## Assigning reward array
slope <- H1 * exp(grid[,4])
for (p in 2:n_position) {
    for (t in 1:(n_dec - 1)) {
        pi <- exp(inflation * (t-1) * step)
        adjust <- discount^(t-1)
        intercept1 <- (slope + H2 * pi) * adjust - (slope * adjust * grid[,4])
        intercept2 <- (slope + (H2 - switch) * pi) * adjust - (slope * adjust * grid[,4])
        if (p > (levels + 1)) {  ## Closed
            ## Close
            reward[, 1, 1, p, t] <- adjust * -maint * pi
            ## Open
            reward[, 1, 3, p, t] <- intercept2
            reward[, 4, 3, p, t] <- slope * adjust
        } else if (p <= (levels + 1)) {  ## Opened
            ## Close
            reward[, 1, 1, p, t] <- -(maint + switch) * pi * adjust
            ## Open
            reward[, 1, 3, p, t] <- intercept1
            reward[, 4, 3, p, t] <- slope * adjust               
        }
    }
}

## Subgrad rep of scrap
scrap <- array(data = 0, dim = c(n_grid, n_dim, n_position))

## Performing fast bellman recursion
time1 <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, disturb_weight, r_index)
time1 <- proc.time() - time1

rm(grid_disturb)
rm(path_sgrid)
rm(reward)
rm(scrap)
rm(disturb)

## Path disturbnaces for duality
n_path <- 500
path_disturb <- array(data = matrix(c(1,     0,     0,     0,
                                      alpha, beta1, beta2, 0,
                                      0,     0,     0,     0,
                                      0,     0,     0,     phi),
                                    ncol = 4, byrow = TRUE),
                      dim = c(n_dim, n_dim, n_path, (n_dec - 1)))
##noise <- rnorm(n_path * (n_dec - 1))
noise <- rnorm(n_path * (n_dec - 1) / 2)
noise <- as.vector(rbind(noise, -noise))  ## anti-thetic disturbances
path_disturb[3, 1,,] <- alpha * noise^2
path_disturb[3, 2,,] <- beta1 * noise^2
path_disturb[3, 3,,] <- beta2 * noise^2
path_disturb[4, 1,,] <- mu * step + alpha * noise * sqrt(step)
path_disturb[4, 2,,] <- beta1 * noise * sqrt(step)
path_disturb[4, 3,,] <- beta2 * noise * sqrt(step)

## Reward function
RewardFunc <- function(state, time) {
    Fpi <- exp(inflation * (time-1) * step)
    Fdiscount <- exp(-(rate + tax_property) * (time-1) * step)
    FH1 <- 5 * step * Fdiscount
    FH2 <- -2.5 * step * Fdiscount * Fpi
    Fmaint <- 0.5 * step *  Fdiscount * Fpi
    Fswitch <- 0.2 * Fdiscount * Fpi
    output <- array(0, dim = c(nrow(state), n_action, n_position))
    output[,1,(2:(levels+1))] <- -(Fmaint + Fswitch)
    output[,3,(2:(levels+1))] <- FH2 + FH1 * exp(state[,4]) 
    output[,1,(levels+2):n_position] <- -Fmaint
    output[,3,(levels+2):n_position] <- FH2 - Fswitch + FH1 * exp(state[,4])       
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(data = 0, dim = c(nrow(state), n_position))
    return(output)
}

## Specifying subsimulation disturbance matrices
n_subsim <- 500
subsim_weight <- rep(1 / n_subsim, n_subsim)
subsim <- array(data = matrix(c(1,     0,     0,     0,
                                        alpha, beta1, beta2, 0,
                                        0,     0,     0,     0,
                                        0,     0,     0,     phi),
                                      ncol = 4, byrow = TRUE),
                        dim = c(n_dim, n_dim, n_subsim, n_path, (n_dec - 1)))
noise <- rnorm(n_subsim * n_path * (n_dec - 1)/ 2)
noise <- as.vector(rbind(noise, -noise))
##noise <- rnorm(n_subsim * n_path * (n_dec - 1))
subsim[3, 1,,,] <- alpha * noise^2
subsim[3, 2,,,] <- beta1 * noise^2
subsim[3, 3,,,] <- beta2 * noise^2
subsim[4, 1,,,] <- mu * step + alpha * noise * sqrt(step)
subsim[4, 2,,,] <- beta1 * noise * sqrt(step)
subsim[4, 3,,,] <- beta2 * noise * sqrt(step)

cat("subsim_done\n")

## Looping
counter <- 1
css_open <- rep(NA, length = nrow(start_index))
bounds_open <- matrix(NA, nrow = nrow(start_index), ncol = 2)
css_closed <- rep(NA, length = nrow(start_index))
bounds_closed <- matrix(NA, nrow = nrow(start_index), ncol = 2)

for (i in 1:nrow(start_index)) {
    temp_start <- start_index[i,]
    path <- PathDisturb(temp_start, path_disturb)
    path_action <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
    ## Primal and dual values
    time2 <- proc.time()
    mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
    duality <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, path_action)
    time2 <- proc.time() - time2
    ## Storing values
    p <- 61
    path_nn <- rflann::FastKDNeighbour(matrix(temp_start, nrow = 1), grid, 1)
    css_open[counter] <- sum(bellman$value[path_nn,,p, 1] * temp_start)
    bounds_open[counter,] <- GetBounds(duality, 0.01, p)
    p <- 121
    css_closed[counter] <- sum(bellman$value[path_nn,,p, 1] * temp_start)    
    bounds_closed[counter,] <- GetBounds(duality, 0.01, p)
    counter <- counter + 1
}

## Results
print(round(cbind(css_open, bounds_open, css_closed, bounds_closed),3))
print(time1)
print(time2)

