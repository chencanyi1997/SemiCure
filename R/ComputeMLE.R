G_gamma <- function(x, gamma) {
  s = exp(-x) * (gamma == 0) + (1 + gamma * x) ^ (-1 / gamma) * (gamma > 0)
  return(s)
}
G_inv <- function(y, gamma) {
  if (gamma == 0) {
    return(-log(y))
  }
  if (gamma > 0) {
    return (y ^ (-gamma) - 1) / gamma
  }
}
#
Inv_F_generic = function(u, t0) {
  a = (1 - u) * (1 - exp(-t0))
  return(-log(exp(-t0) + a))
}

Sfcn = function(gamma, F_hat, beta, Z) {
  return(G_gamma(exp(t(beta) %*% Z) * F_hat, gamma))
}

F.hat = function(alpha) {
  m = length(alpha) + 1
  res = rep(0, m + 1)
  sum = 0
  for (k in 2:(m + 1)) {
    sum = sum + exp(alpha[k - 1])
    res[k] = 1 - exp(-sum)
  }
  res[1] = 0
  res[m + 1] = 1
  return(res)
}

#SEE
#first-order derivative
G_gamma_1d = function(x, gamma) {
  s = -(1 + gamma * x) ^ (-(1 + 1 / gamma)) * (gamma > 0) + (-exp(-x)) * (gamma ==
                                                                            0)
  return(s)
}
#second-order derivative
G_gamma_2d = function(x, gamma) {
  s = (gamma + 1) * (1 + gamma * x) ^ (-(2 + 1 / gamma)) * (gamma > 0) + exp(-x) * (gamma ==
                                                                                      0)
}

Lc_log_generic = function(para, Y, Z, S, Delta, gamma, n, t0) {
  beta = para[1:3]
  alpha = para[-(1:3)]
  Fhat = F.hat(alpha)
  contri = rep(0, n)
  for (i in 1:n) {
    if (Delta[i] == 1) {
      k = which(Y[i] == S)
      contri[i] = log(1 - Sfcn(gamma, Fhat[k + 1], beta, Z[i,]))
    }
    if (Delta[i] == 0) {
      if (Y[i] == t0) {
        k = which(Y[i] == S)
        contri[i] = log(Sfcn(gamma, Fhat[k], beta, Z[i,]))
      }
      else{
        k = which(Y[i] <= S)[1] - 1
        contri[i] = log(Sfcn(gamma, Fhat[k + 1], beta, Z[i,]))
      }
    }
  }
  return(sum(contri))
}

Jobs_generic = function(para, Y, Z, S, Delta, gamma, n, t0) {
  beta1 = para[1:3]
  alpha = para[-(1:3)]
  Fhat = F.hat(alpha)
  J = 0
  for (i in 1:n) {
    # a = exp(sum(t(beta1)*Z[i,])) * F.hat(alpha)

    if (Delta[i] == 1) {
      k = which(Y[i] == S)
      a = exp(sum(t(beta1) * Z[i,])) * Fhat[k + 1]
      a0 = G_gamma(a, gamma)
      a1 = G_gamma_1d(a, gamma)
      a2 = G_gamma_2d(a, gamma)
      Ji = -(a2 * a + a1) * (1 - a0) - a1 ^ 2 * a
      Ji = Ji / (1 - a0) ^ 2 * a * t(t(Z[i,])) %*% Z[i,]
      J = J - Ji
    }
    if (Delta[i] == 0) {
      if (Y[i] == t0) {
        k = which(Y[i] == S)
        a = exp(sum(t(beta1) * Z[i,])) * Fhat[k]
      }
      else{
        k = which(Y[i] <= S)[1] - 1
        a = exp(sum(t(beta1) * Z[i,])) * Fhat[k + 1]
      }
      a0 = G_gamma(a, gamma)
      a1 = G_gamma_1d(a, gamma)
      a2 = G_gamma_2d(a, gamma)
      Ji = (a2 * a + a1) * a0 - a1 ^ 2 * a
      Ji = Ji / a0 ^ 2 * a * t(t(Z[i,])) %*% Z[i,]
      J = J - Ji
    }
  }
  return(J)
}

#' Compute_MLE
#'
#' @param n number of observations
#' @param gamma rate
#' @param beta true beta
#' @param t0 censor time
#'
#' @return estimated beta
#' @export
#'
#' @examples
#' library(doSNOW)
#' rep <- 10
#' n <- 200
#' beta <- c(-1 / 2, 1,-1 / 2)
#' t0 <- 4
#' gamma <- 0.75
#' cl <- makeCluster(10)
#' registerDoSNOW(cl)
#' result = foreach(i = 1:rep, .packages = "semicure") %dopar% {
#'    .libPaths(c(.libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.6"))
#'    n <- 200
#'    beta <- c(-1 / 2, 1,-1 / 2)
#'    t0 <- 4
#'    gamma <- 0.75
#'    semicure::Compute_MLE(n, gamma, beta, t0)
#' }
#' stopCluster(cl)
#'
#'
#' #statistic
#' estimate = matrix(0, ncol = rep, nrow = 2)
#' SEE = matrix(0, ncol = rep, nrow = 2)
#' for (i in 1:rep) {
#'   estimate[, i] = (result[[i]][[1]])[-1]
#'   SEE[, i] = (result[[i]][[2]])[-1]
#' }
#' rowMeans(estimate) # estimate
#' sqrt(cov(t(estimate))) # SE
#' rowMeans(SEE) #SEE
#'
#' #CP ?????????????
#' c = qnorm(1 - 0.025)
#' cnt = c(0, 0)
#' for (i in 1:rep) {
#'   Test = abs(estimate[, i] - beta[-1]) <= c * SEE[, i]
#'   cnt = cnt + Test
#' }
#' CP = cnt / rep
#' CP
Compute_MLE = function(n, gamma, beta, t0) {
  Inv_F <- pryr::partial(Inv_F_generic, t0 = t0)
  # n = 400
  # gamma = 0
  Z1 = runif(n)
  Z2 = rbinom(n, 1, 0.5)
  Z = cbind(1, Z1, Z2)

  #generate T
  T <- rep(NA, n)
  for (i in 1:n) {
    u = runif(1)
    if (u > (1 - G_gamma(exp(t(beta) %*% Z[i,]), gamma))) {
      T[i] = Inf
    }
    if (u <= (1 - G_gamma(exp(t(beta) %*% Z[i,]), gamma))) {
      T[i] = Inv_F(exp(-t(beta) %*% Z[i,]) * G_inv(1 - u, gamma))
    }
  }
  # length(which(T>1e4))/n # cure rate

  Y = rexp(n, 0.5)
  Y = pmin(t0, Y)
  Delta = (T <= Y) * 1

  S = Y[which(Delta == 1)]
  S = unique(S)
  S = sort(S)
  m = length(S)

  ### Initial value for alpha(jump size = 1/m at S1,...Sm)
  alpha_init = rep(0, m - 1)
  alpha_init[1] = log(-log(1 - 1 / m))
  sum = 0
  for (j in 2:(m - 1)) {
    sum = sum + exp(alpha_init[j - 1])
    alpha_init[j] = log(-log(1 - exp(sum) / m))
  }

  Lc_log <-
    pryr::partial(
      Lc_log_generic,
      Y = Y,
      Z = Z,
      S = S,
      Delta = Delta,
      gamma = gamma,
      n = n,
      t0 = t0
    )
  Jobs <-
    pryr::partial(
      Jobs_generic,
      Y = Y,
      Z = Z,
      S = S,
      Delta = Delta,
      gamma = gamma,
      n = n,
      t0 = t0
    )

  res = optim(
    par = c(0, 0, 0, alpha_init),
    fn = Lc_log,
    method = "BFGS",
    control = list(fnscale = -1)
  )

  para = res$par
  SEE = Jobs(para)
  SEE = sqrt(solve(SEE))
  return(list(para[1:3], diag(SEE)))
}

#run time
# microbenchmark(Compute_MLE(200,0), times = 1)


# #parallel

#run time
# microbenchmark(Compute_MLE(200,0), times = 1)

# result <- Compute_MLE(200,0)
# profvis(Compute_MLE(200,0))
# microbenchmark::microbenchmark(Compute_MLE(200,0),times = 1)
