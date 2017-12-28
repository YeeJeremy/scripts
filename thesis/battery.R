## Battery control

rm(list = ls())
gc()
library(rcss)

## Grid
n_grid <- 501  ## number of grid points
grid_start <- -15  ## lowest state
grid_end <- 15  ## highest state
grid <- cbind(rep(1, n_grid), seq(grid_start, grid_end, length = n_grid))

## Battery specification
step <- 5  ## step between battery levels
pos_start <- 0
pos_end <- 100
position <- seq(pos_start, pos_end, by = step)  ## battery levels
n_pos <- length(position)

## Standard deviation for the consumer demand
std <- 10

## Actions
n_action <- 11  ## number of safety margins
safety_start <- 0
safety_end <- 50
safety <- seq(safety_start, safety_end, length = n_action)  ## safety margins

## Control array
control <- array(data = 0, dim = c(n_pos, n_action, n_pos))
for (p in 1:n_pos) {
    for (a in 1:n_action) {
        temp <- position[p] + safety[a]  ## center of normal distribution
        control[p,a,1] <- pnorm(pos_start + step/2, temp, std)
        control[p,a,n_pos] <- 1 - pnorm(pos_end - step/2, temp, std)
        for (pp in 2:(n_pos-1)) {
            control[p,a,pp] <- pnorm(position[pp] + step/2, temp, std) -
                pnorm(position[pp] - step/2, temp, std)
        }
    }
}
    
## Functions to calculate expected excess and shortage energy demand
erf <- function(x){  ## error function
    return(2 * pnorm(x * sqrt(2)) - 1)
}
Excess <- function(pos, act) {
    temp1 <- pos_end + step/2
    temp2 <- pos + act
    result <- std/sqrt(2*pi) * exp(-(temp1-temp2)^2/(2*std^2)) +
        (temp2 - pos_end)/2 * (1 - erf(1/sqrt(2*std^2) * (temp1 - temp2)))
    return(result)
}
Shortage <- function(pos, act) {
    temp1 <- pos_start - step/2
    temp2 <- pos + act
    result <- std/sqrt(2*pi) * exp(-(temp1-temp2)^2/(2*std^2)) +
        (pos_start - temp2)/2 * (erf(1/sqrt(2*std^2) * (temp1 - temp2)) + 1)
    return(result)
}

## Expected excess and shortage energy demand
excess <- matrix(data = NA, nrow = n_pos, ncol = n_action)
shortage <- matrix(data = NA, nrow = n_pos, ncol = n_action)
for (p in 1:n_pos) {
    for (a in 1:n_action) {
        excess[p,a] <- Excess(position[p], safety[a])
        shortage[p,a] <- Shortage(position[p], safety[a])
    }
}

## Subgradient representation of reward functions
n_dec <- 48 * 7 ## number of decision epochs
u_t <- 10 + cos((0:(n_dec-1)) * 2*pi/48 + 3*pi/2)
v_t <- 1 + (sin((0:(n_dec-1)) * 2*pi/48 + 3*pi/2))/2
buy <- 20  ## price to buy from grid
sell <- 0  ## price to sell to grid
reward <- array(0, dim = c(n_grid, 2, n_action, n_pos, n_dec - 1))
for (p in 1:n_pos) {
    for (a in 1:n_action) {
        for (t in 1:(n_dec-1)) {
            reward[,1,a,p,t] <- -safety[a] * u_t[t] - shortage[p, a] * buy + excess[p, a] * sell
            reward[,2,a,p,t] <- -safety[a] * v_t[t]
        }
    }
}
scrap <- array(0, dim = c(n_grid, 2, n_pos))
for (p in 1:n_pos) {
    scrap[,1,p] <- position[p] * u_t[n_dec]
    scrap[,2,p] <- position[p] * v_t[n_dec]
}

## Parameters for AR(1) process (Z_t)
mu <- 0
sigma <- 0.5
phi <- 0.9

## Disturbance sampling
n_disturb <- 1000  ## size of sampling
disturb_weight <- rep(1/n_disturb, n_disturb)  ## probability weights
disturb <- array(matrix(c(1, 0, 0, phi), ncol = 2, byrow = TRUE), dim = c(2, 2, n_disturb))
CondExpected <- function(a, b){
    return(1/sqrt(2 * pi) * (exp(-a^2/2)- exp(-b^2/2)))
}
part <- qnorm(seq(0, 1, length = n_disturb + 1))
for (i in 1:n_disturb) {
    disturb[2,1,i] <- mu + sigma * (CondExpected(part[i], part[i+1]) / (pnorm(part[i+1]) - pnorm(part[i])))
}
r_index <- matrix(c(2, 1), ncol = 2)  ## randomness index

## Fast bellman recursion
time1 <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, disturb_weight, r_index)
time1 <- proc.time() - time1

## Solution Diagnostics
## Exact reward function
RewardFunc <- function(state, time) {
    output <- array(0, dim = c(nrow(state), n_action, n_pos))
    for (p in 1:n_pos) {
        for (a in 1:n_action) {
            output[,a,p] <- -safety[a] * (u_t[time] + v_t[time] * state[,2]) - shortage[p,a] * buy + excess[p,a] * sell
        }
    }
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(0, dim = c(nrow(state), n_pos))
    for (p in 1:n_pos) {
        output[,p] <- position[p] * (u_t[n_dec] + v_t[n_dec] * state[,2])
    }
    return(output)
}

## Generate sample path disturbances
set.seed(12345)
n_path <- 100
path_disturb <- array(matrix(c(1, 0, 0, phi), ncol = 2, byrow = TRUE),
                      dim = c(2, 2, n_path, n_dec - 1))
rand <- rnorm(n_path * (n_dec - 1) / 2)
rand <- as.vector(rbind(rand, -rand))
path_disturb[2,1,,] <- mu + sigma * rand
start <- c(1, 0)  ## z_0
path <- PathDisturb(start, path_disturb)
policy <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)

## Specifying subsimulation disturbance matrices
n_subsim <- 100
subsim_weight <- rep(1/n_subsim, n_subsim)
subsim <- array(matrix(c(1, 0, 0, phi), ncol = 2, byrow = TRUE), dim = c(2, 2, n_subsim, n_path, n_dec - 1))
rand <- rnorm(n_subsim * n_path * (n_dec - 1)/2)
rand <- as.vector(rbind(rand, -rand))
subsim[2,1,,,] <- mu + sigma * rand

## Compute primal and dual values
time2 <- proc.time()
mart <- FastAddDual(path, subsim, subsim_weight, grid, bellman$value, ScrapFunc)
bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, policy)
time2 <- proc.time() - time2

## Storing the results
results <- matrix(data = NA, nrow = n_pos, ncol = 3)
alpha <- 0.01
for (p in 1:n_pos) {
    results[p,1] <- sum(bellman$value[251,,p,1] * grid[251,])
    results[p,2:3] <- GetBounds(bounds, alpha, p)
}

## Results
print(round(results,3))
print(time1)
print(time2)
