## Testing the choice of sampling
rm(list = ls())
gc()
library(rcss)

## Grid of evenly spaced points
n_grid <- 301 
grid <- as.matrix(cbind(rep(1, n_grid), seq(30, 60, length = n_grid)))

## Parameters
rate <- 0.06
step <- 0.02
vol <- 0.2
n_dec <- 51

## Transition
control <- matrix(c(c(1, 1), c(2, 1)), nrow = 2, byrow = TRUE)

## Subgrad rep of reward
strike <- 40
in_money <- grid[,2] <= strike
reward <- array(0, dim = c(n_grid, 2, 2, 2, n_dec - 1))       
reward[in_money, 1, 2, 2,] <- strike
reward[in_money, 2, 2, 2,] <- -1
for (tt in 1:(n_dec - 1)){
  reward[,,,,tt] <- exp(-rate * step * (tt - 1)) * reward[,,,,tt] 
}

## Subgrad rep of scrap
scrap <- array(data = 0, dim = c(n_grid, 2, 2))
scrap[in_money, 1, 2] <- strike
scrap[in_money, 2, 2] <- -1
scrap <- exp(-rate * step * (n_dec - 1)) * scrap

## Disturbance parameters
r_index <- matrix(c(2, 2), ncol = 2)
u <- (rate - 0.5 * vol^2) * step
sigma <- vol * sqrt(step)
condExpected <- function(a, b){  ##[a,b]
    aa <- (log(a) - (u+sigma^2))/sigma
    bb <- (log(b) - (u+sigma^2))/sigma
    vv <- exp(u + sigma^2/2) * (pnorm(bb) - pnorm(aa))
    return(vv)
}

## Container
valueLA <- rep(NA, 10)
valueMC <- valueLA
ddIndex <- c(seq(1000, 10000, by = 1000))
ind <- which(grid[,2] == 36)
set.seed(12345)

for (dd in 1:length(ddIndex)) {
    ## Disturbance sampling
    n_disturb <- ddIndex[dd]
    weight <- rep(1/n_disturb, n_disturb)
    disturb <- array(0, dim = c(2, 2, n_disturb))
    disturb[1, 1,] <- 1
    ## Local averages
    quantile <- rep(NA, n_disturb)
    part <- qlnorm(seq(0, 1, length = n_disturb + 1), u, sigma)
    for (i in 1:n_disturb) {
        quantile[i] <- condExpected(part[i], part[i+1]) / (plnorm(part[i+1], u, sigma) - plnorm(part[i], u, sigma))
    }
    disturb[2, 2,] <- quantile
    ## Performing fast bellman recursion
    bellmanLA <- FastBellman(grid, reward, scrap, control, disturb, weight, r_index)
    valueLA[dd] <- sum(bellmanLA$value[ind, , 2, 1] * grid[ind,])
    ## Use random samples    
    rand1 <- rnorm(n_disturb/2, 0, 1)
    rand1 <- c(rand1, -rand1)
    rand1 <-  exp((rate - 0.5 * vol^2) * step + vol * sqrt(step) * rand1)
    disturb[2, 2,] <- rand1
    ## Performing fast bellman recursion
    bellmanMC <- FastBellman(grid, reward, scrap, control, disturb, weight, r_index)
    valueMC[dd] <- sum(bellmanMC$value[ind, , 2, 1] * grid[ind,])
}

## Plot
setEPS()
postscript("BputC.eps")
plot(ddIndex, valueMC, type = "l", xlab = "Sample Size", ylab = "Value Estimate")
dev.off()
dev.new()
setEPS()
postscript("BputLA.eps")
plot(ddIndex, valueLA, type ="l", xlab = "Sample Size", ylab ="Value Estimate")
dev.off()
