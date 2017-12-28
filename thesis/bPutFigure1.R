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
disturb[1, 1,] <- 1
part <- qlnorm(seq(0, 1, length = 1000 + 1), u, sigma)
for (i in 1:1000) {
    disturb[2, 2, i] <- condExpected(part[i], part[i+1]) / (plnorm(part[i+1], u, sigma) - plnorm(part[i], u, sigma))
}

## Looping
ggIndex <- 1:10
times <- matrix(NA, nrow = 10, ncol = 3)
diffs <- matrix(NA, nrow = 10, ncol = 2)

for (gg in 1:length(ggIndex)) {
    ## Grid
    n_grid <- ggIndex[gg] * 30 + 1
    grid <- as.matrix(cbind(rep(1, n_grid), seq(30, 60, length = n_grid)))
    ## Subgradient representation of reward
    in_money <- grid[,2] <= strike
    reward <- array(0, dim = c(n_grid, 2, 2, 2, n_dec))       
    reward[in_money, 1, 2, 2,] <- strike
    reward[in_money, 2, 2, 2,] <- -1
    for (tt in 1:n_dec){
        reward[,,,,tt] <- exp(-rate * step * (tt - 1)) * reward[,,,,tt]
    }
    ## Subgrad representation of scrap
    scrap <- array(data = 0, dim = c(n_grid, 2, 2))
    scrap[in_money, 1, 2] <- strike
    scrap[in_money, 2, 2] <- -1
    scrap <- exp(-rate * step * (n_dec - 1)) * scrap
    ## Different
    time <- proc.time()
    bellman <- Bellman(grid, reward, scrap, control, disturb, weight)
    time <- proc.time() - time
    time1 <- proc.time()
    bellman1 <- AcceleratedBellman(grid, reward, scrap, control, disturb, weight, 2)
    time1 <- proc.time() - time1
    time2 <- proc.time()
    bellman2 <- FastBellman(grid, reward, scrap, control, disturb, weight, r_index)
    time2 <- proc.time() - time2
    ## Value function approximations
    value <- rowSums(bellman$value[,,2,1] * grid)
    value1 <- rowSums(bellman1$value[,,2,1] * grid)
    value2 <- rowSums(bellman2$value[,,2,1] * grid)
    ## Get absolute difference
    diffs[gg,1] <- max(abs(value1 - value))
    diffs[gg,2] <- max(abs(value2 - value))
    ## Store values
    times[gg, 1] <- time[1]
    times[gg, 2] <- time1[1]
    times[gg, 3] <- time2[1]    
}

setEPS()
postscript("Bput_nn_time.eps")
x <- seq(31, 301, 30)
matplot(x, times, type = "l", main = "CPU Time", xlab = "Grid Points", ylab = "")
dev.off()

setEPS()
postscript("Bput_nn_diff.eps")
matplot(x, diffs, type = "l", main = "Max Difference", xlab = "Grid Points", ylab ="")
dev.off()
