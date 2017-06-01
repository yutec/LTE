heidelberger <- function(dat1,thin=1){
  ndim = 8
  epsilon = 0.05
  pval = 0.10
  #colnames(stat) = c("sdBeta1","sdBeta2","sdAlpha","beta0","beta1","beta2","beta3","alpha","Init","k","Time")

  if (any(is.na(dat1))==TRUE){
    stat = c(0,-200)
  }
  else if (sd(dat1[dim(dat1)[1]-(99:0),1])==0){
    stat = c(0,-300)
  }
  else {
    print(paste("Heidelberger test for the latest draws of ",dim(dat1)[1]))
    x = seq(1, dim(dat1)[1], by=thin)
    h = heidel.diag(dat1[x,], eps=epsilon, pvalue=pval)
    rownames(h) = c("sdBeta1","sdBeta2","sdAlpha","beta0","beta1","beta2","beta3","alpha")
    print(h)

    flag = (min(h[,4])!=1 | max(is.na(h[,4]))==1)
    stat = c(max(h[,2]),-100*flag)
    stat[is.na(stat)] = 0
  }
  return(stat)

  print("Heidelberger complete.")
}

heidel1 <- function(x, pvalue=0.10, thin=1) {
  x <- as.mcmc(as.matrix(x))
  x <- seq(1, dim(x)[1], by=thin)
  HW.mat0 <- matrix(0, ncol = 3, nrow = nvar(x))
  dimnames(HW.mat0) <- list(varnames(x), c("stest", "start", "pvalue"))
  HW.mat <- HW.mat0
  for (j in 1:nvar(x)) {
    start.vec <- seq(from = start(x), to = end(x)/2, by = niter(x)/10)
    Y <- x[, j, drop = TRUE]
    n1 <- length(Y)
    S0 <- spectrum0.ar(window(Y, start = end(Y)/2))$spec
    converged <- FALSE
    for (i in seq(along = start.vec)) {
      Y <- window(Y, start = start.vec[i])
      n <- niter(Y)
      ybar <- mean(Y)
      B <- cumsum(Y) - ybar * (1:n)
      Bsq <- (B * B)/(n * S0)
      I <- sum(Bsq)/n
      if (converged <- !is.na(I) && pcramer(I) < 1 - pvalue)
        break
    }
    if (!converged || is.na(I)) {
      nstart <- NA
    } else {
      nstart <- start(Y)
    }
    HW.mat[j, ] <- c(converged, nstart, 1 - pcramer(I))
  }
  return(HW.mat)
}

heidel2 <- function(x, eps=0.05, thin=1) {
  x <- as.mcmc(as.matrix(x))
  x <- seq(1, dim(x)[1], by=thin)
  HW.mat0 <- matrix(0, ncol = 3, nrow = nvar(x))
  dimnames(HW.mat0) <- list(varnames(x), c("htest", "mean", "halfwidth"))
  HW.mat <- HW.mat0

  for (j in 1:nvar(x)){
    Y <- x[, j, drop = TRUE]
    #n <- niter(Y)
    n <- length(Y)
    ybar <- mean(Y)
    S0ci <- spectrum0.ar(Y)$spec
    halfwidth <- 1.96 * sqrt(S0ci/n)
    passed.hw <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
    HW.mat[j, ] <- c(passed.hw, ybar, halfwidth)
  }
  return(HW.mat)
}

htest <- function(dat1){
  ndim = 8
  epsilon = 0.05
  pval = 0.10
  h = heidel2(dat1, eps=epsilon)
  rownames(h) = c("sdBeta1","sdBeta2","sdAlpha","beta0","beta1","beta2","beta3","alpha")
  print(h)
  flag = (min(h[,1])!=1 | max(is.na(h[,1]))==1)
  stat = c(max(h[,2]),-100*flag)
  stat[is.na(stat)] = 0
  return(stat)
}
