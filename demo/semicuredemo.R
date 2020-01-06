library(doSNOW)
rep <- 10
n <- 200
beta <- c(-1 / 2, 1,-1 / 2)
t0 <- 4
gamma <- 0.75
cl <- makeCluster(4)
registerDoSNOW(cl)
result = foreach(i = 1:rep, .packages = "semicure") %dopar% {
  # brary("semicure", lib.loc = "/home/chencanyi2018/R/x86_64-pc-linux-gnu-library/3.6")
  n <- 200
  beta <- c(-1 / 2, 1,-1 / 2)
  t0 <- 4
  gamma <- 0.75
  s <- Compute_MLE(n = n, gamma = gamma, beta = beta, t0 = t0)
  return(s)
}
stopCluster(cl)


#statistic
estimate = matrix(0, ncol = rep, nrow = 2)
SEE = matrix(0, ncol = rep, nrow = 2)
for (i in 1:rep) {
  estimate[, i] = (result[[i]][[1]])[-1]
  SEE[, i] = (result[[i]][[2]])[-1]
}
rowMeans(estimate) # estimate
sqrt(cov(t(estimate))) # SE
rowMeans(SEE) #SEE

#CP ?????????????
c = qnorm(1 - 0.025)
cnt = c(0, 0)
for (i in 1:rep) {
  Test = abs(estimate[, i] - beta[-1]) <= c * SEE[, i]
  cnt = cnt + Test
}
CP = cnt / rep
CP
