## Point and figure charting table of results. Second table.

rm(list = ls())
gc()
library(rcss)
set.seed(12345)                   

## Parameters
n_dec <- 101                        ## Number of decision epochs
step <- 1                          ## Price change step
cost <- 0.5                        ## Cost to change positions 

## Transition kernels
p1 <- 0.8                          ## Prob(state 1 -> 1)
p2 <- 0.8                          ## Prob(state 2 -> 2)
q1 <- 0.9                          ## Prob(signal = 1 | state = 1)
q2 <- 0.9                          ## Prob(signal = 2 | state = 2)

## Transition matrix
n_state <- 2                       ## Number of hidden state estimates
gamma <- matrix(c(p1, 1-p2, 1-p1, p2), nrow = n_state, ncol = n_state)

## Reference measure
n_disturb <- 2                     ## Number of observation types
ref <- c(0.5, 0.5)                 ## Reference measure

## Diagnonal entries representing conditional probabilities for observations
diag_obs1 <- c(q1, 1-q1)           ## Observation 1
diag_obs2 <- c(1-q2, q2)           ## Observation 2

## Control matrix
n_pos <- 3                         ## Number of positions
n_action <- 3                      ## Number of actions
control <- matrix(c(1,1,1,2,2,2,3,3,3), nrow = n_pos, ncol = n_action)

## Disturbance sampling
disturb <- array(0, dim = c(n_state + 1, n_state + 1, n_disturb))
disturb[1,1,] <- 1
disturb[-1,-1,1] <- t(gamma) %*% diag(as.vector(diag_obs1/ref[1]))
disturb[-1,-1,2] <- t(gamma) %*% diag(as.vector(diag_obs2/ref[2]))
r_index <- matrix(c(2,2,2,3,3,2,3,3), ncol = 2 , byrow =TRUE)

## Equally spaced grid
temp_n_grid <- 6  ## Size of the seed grid
seed_grid <- cbind(rep(1, temp_n_grid), seq(0, 1, length = temp_n_grid))
seed_grid <- cbind(seed_grid, 1 - seed_grid[,2])
grid <- seed_grid
counter1 <- 0
counter2 <- temp_n_grid
for (t in 1:(9 - 1)) {
    count <- 2^(t-1)  ## number of grids to use
    for (d in 1:count) {
        index <- counter1 + (d - 1) * temp_n_grid + 1 ## starting row
        temp1 <- grid[(index):(index + temp_n_grid - 1),] %*% t(disturb[,,1])
        temp2 <- grid[(index):(index + temp_n_grid - 1),] %*% t(disturb[,,2])
        grid <- rbind(grid, temp1, temp2)
    }
    counter1 <- counter2  ## previous number in grid
    counter2 <- counter1 + count * temp_n_grid * 2 ## number in grid
}
n_grid <- nrow(grid)

library(lattice)
x <- seq(0,40,by=1)
y <- seq(0,40,by=1)
grid1 <- expand.grid(x=x, y=y)
n_grid1 <- nrow(grid1)
grid1 <- as.matrix(cbind(rep(1,n_grid1),grid1))
grid <- as.matrix(rbind(grid, grid1))
n_grid <- nrow(grid)

## Subgradient representation of reward function
reward <- array(0, dim = c(n_grid, n_state + 1, n_action, n_pos, n_dec - 1))
for (p in 1:n_pos) {
    for (a in 1:n_action) {
        for (tt in 1:(n_dec - 1)) {
            ## The expected return from portfolio
            reward[,2,a,p,tt] <- (control[p,a] - 2) * step * (1 - 2 * q1)
            reward[,3,a,p,tt] <- (control[p,a] - 2) * step * (2 * q2 - 1)
            ## Transaction cost
            reward[,1,a,p,tt] <- -cost * abs(p - control[p,a])
        }
    }
}
## Sudbgradient representation of scrap
scrap <- array(0, dim = c(n_grid, n_state + 1, n_pos))
for (p in 1:n_pos) {
    scrap[,1,p] <- -cost * abs(p - control[p,2])
}

## Bellman recursion
time1 <- proc.time()
bellman <- FastBellman(grid, reward, scrap, control, disturb, ref, r_index)
time1 <- proc.time() - time1

## Solution diagnostic
n_path <- 50000

## Simulate path disturbances
path_disturb <- array(0, dim=c(3, 3, n_path, (n_dec - 1)))
for (i in 1:(n_dec - 1)) {
    for (j in 1:(n_path / 2)) {
        index <- sample(c(1, 2), 1)
        path_disturb[,, j, i] <- disturb[,,index]
        path_disturb[,, n_path/2 + j, i] <- disturb[,,index%%2 + 1]
    }
}

## Reward function
RewardFunc <- function(state, time) {
    output <- array(0, dim = c(nrow(state), n_action, n_pos))
    change <- step * ((1 - 2 * q1) * state[,2] + (2 * q2 -1) * state[,3])
    for (p in 1:n_pos) {
        for (a in 1:n_action) {
            output[,a,p] <- (a - 2) * change - cost * abs(p - control[p,a])
        }
    }
    return(output)
}
## Scrap function
ScrapFunc <- function(state) {
    output <- array(0, dim = c(nrow(state), n_pos))
    for (p in 1:n_pos) {
        output[,p] <- -cost * abs(p - control[p,2])  
    }
    return(output)
}

## Subsimulation disturbances
subsim <- array(0, dim = c(3, 3, n_disturb, n_path, n_dec - 1))
for (i in 1:(n_dec-1)) {
    for (j in 1:n_path) {
        for (k in 1:n_disturb) {
            subsim[,,k,j,i] <- disturb[,,k]
        }
    } 
}

## Run the diagnostic check and testing for various starting points
start_index <- seq(0, 1, by = 0.2)
css <- rep(NA, length(start_index) * 3)
test <- css
test_se <- css
ci <- matrix(NA, nrow = (length(start_index) * 3), ncol = 2)
count <- 0
for (ss in start_index) {
    start <- c(1, ss, 1 - ss)
    path <- PathDisturb(start, path_disturb)
    path_action <- FastPathPolicy(path, grid, control, RewardFunc, bellman$expected)
    time2 <- proc.time()
    mart <- FastAddDual(path, subsim, ref, grid, bellman$value, ScrapFunc)
    bounds <- AddDualBounds(path, control, RewardFunc, ScrapFunc, mart, path_action)
    time2 <- proc.time() - time2
    path_nn <- rflann::FastKDNeighbour(matrix(start, nrow = 1), grid, 1)
    for (i in 1:3) { 
        css[count + i] <- sum(bellman$value[path_nn, ,i, 1] * start)
        temp_test <- TestPolicy(i, path, control, RewardFunc, ScrapFunc, path_action)
        test[count + i] <- mean(temp_test)
        test_se[count + i] <- sd(temp_test)/sqrt(n_path)
        ci[count + i, ] <- GetBounds(bounds, 0.01, i)
    }
    count <- count + 3
}

## Displaying results
table <- cbind(css, test, test_se,ci)
print(round(table, 5))
print(time1)
print(time2)

## Plot the grid used
plot(grid[,2:3], xlab ="Probability of State 1", ylab = "Probability of State 2")

## Plot the value functions
plotNGrid <- 101
plotGrid <- cbind(rep(1, plotNGrid), seq(0, 1, length = plotNGrid))
plotGrid <- cbind(plotGrid, 1 - plotGrid[,2])
valueF1 <- Optimal(plotGrid, bellman$value[,,1,1])
valueF2 <- Optimal(plotGrid, bellman$value[,,2,1])
valueF3 <- Optimal(plotGrid, bellman$value[,,3,1])
value1 <- rowSums(valueF1 * plotGrid)
value2 <- rowSums(valueF2 * plotGrid)
value3 <- rowSums(valueF3 * plotGrid)
plot(plotGrid[,2], value1, lty=1, type = "l", xlab = "Probability of State 1", ylab = "Value")
lines(plotGrid[,2], value2, lty =2, type = "l", col ="red")
lines(plotGrid[,2], value3, lty =3, type = "l", col = "blue")
