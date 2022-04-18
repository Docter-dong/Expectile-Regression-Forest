########################################################
# Description:
# 1.Purpose: functions for ernn with L1 and (L1+L2)
# 4.Author: Ming J
# 5.Founded: Jan. 3, 2018.
# 6.Revised: Jan. 3, 2018.
########################################################

#######
# 1. ernn.fit
# attention: alpha = 1 means "L2"; alpha = 0 means "L1"
#######
ernn.fit <- function (x, y, n.hidden, tau = 0.5, n.ensemble = 1, iter.max = 5000,  
                      n.trials = 3, bag = FALSE, lower = -Inf,  
                      Th = sigmoid, Th.prime = sigmoid.prime, penalty = 0, trace = TRUE, 
                      eps = 10^-4, alpha = 1, ...) 
{
  if (!is.matrix(x)) 
    stop("\"x\" must be a matrix")
  if (!is.matrix(y)) 
    stop("\"y\" must be a matrix")
  if (ncol(y) != 1) 
    stop("\"y\" must be univariate")
  if ((tau > 1) | (tau < 0)) 
    stop("invalid \"tau\"")
  if (identical(Th, linear))   
    n.hidden <- 1
  x <- scale(x) 
  x.center <- attr(x, "scaled:center")  
  x.scale <- attr(x, "scaled:scale")   
  y <- scale(y)
  y.center <- attr(y, "scaled:center")  
  y.scale <- attr(y, "scaled:scale")
  lower.scaled <- (lower - y.center)/y.scale  
  weights <- list()
  if (trace) 
    cat("tau=", tau, "\n", sep = "")
  for (i in 1:n.ensemble) {
    if (trace) 
      cat(i, "/", n.ensemble, "\n", sep = "")
    w.tmp <- NA
    class(w.tmp) <- "try-error"  
    
    n.errors <- 0
    while (inherits(w.tmp, "try-error")) {  
      w.tmp <- try(ernn.nlm(x, y, n.hidden, tau, iter.max, 
                            n.trials, bag, lower.scaled, Th, Th.prime, 
                            penalty, trace, eps, alpha, ...), silent = TRUE)  
      n.errors <- n.errors + 1   
      if (n.errors > 5) 
        stop("nlm optimization failed")
    }
    weights[[i]] <- w.tmp 
  }
  if (trace) 
    cat("\n")
  parms <- list(weights = weights, lower = lower, 
                tau = tau, Th = Th, x.center = x.center, x.scale = x.scale, 
                y.center = y.center, y.scale = y.scale)
  parms
}


######
# 2. ernn.nlm
######
ernn.nlm <- function (x, y, n.hidden, tau, iter.max, n.trials, bag, lower, 
                      Th, Th.prime, penalty, trace, eps, alpha, ...) 
{
  cases <- 1:nrow(x)
  if (bag) 
    cases <- sample(nrow(x), replace = TRUE)
  x <- x[cases, , drop = FALSE]  
  y <- y[cases, , drop = FALSE]
  if (length(lower) > 1)
    lower <- lower[cases] 
  
  cost.best <- Inf
  for (i in 1:n.trials) {
    weights <- ernn.initialize(x, y, n.hidden)
    print(weights)
    if (any(lower != -Inf)) {  
      
      fit <- nlm(ernn.cost, weights, iterlim = iter.max,
                 x = x, y = y, n.hidden = n.hidden, tau = tau,
                 eps = eps, lower = -Inf, Th = Th, Th.prime = Th.prime, 
                 penalty = penalty, alpha = alpha, check.analyticals = FALSE, 
                 ...)
      weights <- fit$estimate
      
    }
    else {  
      fit <- nlm(ernn.cost, weights, iterlim = iter.max, 
                 x = x, y = y, n.hidden = n.hidden, tau = tau, 
                 eps = eps, lower = lower, Th = Th, Th.prime = Th.prime, 
                 penalty = penalty, alpha = alpha, check.analyticals = FALSE, 
                 ...)
      weights <- fit$estimate
    }
    cost <- fit$minimum
    if (trace) 
      cat(i, cost, "\n")
    if (cost < cost.best) {
      cost.best <- cost
      weights.best <- fit$estimate
    }
  }
  if (trace) 
    cat("*", cost.best, "\n")
  weights.best <- ernn.reshape(x, y, weights.best, n.hidden)
  weights.best
}


######
# 3.ernn.cost
######
ernn.cost <- function (weights, x, y, n.hidden, tau, eps, lower, Th, Th.prime, 
                       penalty, alpha) 
{
  penalty2 <- ifelse(identical(Th, linear), penalty, 0)
  w <- ernn.reshape(x, y, weights, n.hidden)
  W1 <- w$W1
  rW1 <- nrow(W1)
  cW1 <- ncol(W1)
  W2 <- w$W2
  rW2 <- nrow(W2)
  cW2 <- ncol(W2)
  x <- cbind(x, 1)
  h1 <- x %*% W1
  y1 <- Th(h1)
  aug.y1 <- cbind(y1, 1)
  h2 <- aug.y1 %*% W2
  
  y2 <- fexp(h2, lower)
  
  E <- y - y2      
  delta2 <- fexp.prime(h2, lower) * tilted.expe.prime(E,tau)
  # browser()
  gradient.W2 <- -(t(aug.y1) %*% delta2)/length(E) + penalty2*(alpha*2*rbind(W2[1:(rW2 - 1), , drop = FALSE], 0) + 
                                                                 (1 - alpha)*rbind(huber.prime(W2[1:(rW2 - 1), , drop = FALSE], eps), 0) )/(length(W2) - cW2)
  
  E1 <- delta2 %*% t(W2[1:(rW2 - 1), , drop = FALSE])
  
  delta1 <- Th.prime(h1) * E1
  
  gradient.W1 <- -(t(x) %*% delta1)/length(E) + penalty*(alpha*2*rbind(W1[1:(rW1 - 1), , drop = FALSE], 0) + 
                                                           (1-alpha)*rbind(huber.prime(W1[1:(rW1 - 1), , drop = FALSE], eps), 0))/(length(W1) - cW1)
  
  cost <- sum(tilted.expe(E, tau))/length(E) + penalty*
    (alpha*sum(W1[1:(rW1 - 1), , drop = FALSE]^2) + (1-alpha)*sum(huber(W1[1:(rW1 - 1), , drop = FALSE], eps)))/(length(W1)-cW1) + 
    penalty2*(alpha*sum(W2[1:(rW2 - 1), , drop = FALSE]^2) + (1-alpha)*sum(huber(W2[1:(rW2 - 1), , drop = FALSE], eps)))/(length(W2) - cW2)
  
  gradient <- c(gradient.W1, gradient.W2)
  attr(cost, "gradient") <- gradient
  cost
}

######
# 4.ernn.reshape
######
ernn.reshape <- function (x, y, weights, n.hidden)
{
  N11 <- ncol(x) + 1
  N12 <- n.hidden
  N1 <- N11 * N12
  W1 <- weights[1:N1]
  W1 <- matrix(W1, N11, N12)
  N21 <- n.hidden + 1
  N22 <- ncol(y)
  N2 <- N1 + N21 * N22
  W2 <- weights[(N1 + 1):N2]
  W2 <- matrix(W2, N21, N22)
  list(W1 = W1, W2 = W2)
}

######
# 5. Th = sigmoid; Th.prime = sigmoid.prime
######
sigmoid <- function(x)
{
  ans <- 1/(1+exp(-x))
  ans
}

sigmoid.prime <- function (x)
{
  ans <- sigmoid(x)*(1 - sigmoid(x))
  ans
}

linear <- function (x) 
{
  x
}
######
# 6.fexp
######
fexp <- function (x, lower)
{
  if (length(lower) > 1) { 
    mapply(fexp, x, lower) 
  }
  else {
    if (lower == -Inf) {
      return(x)
    }
    else {
      stop("\"lower\" should be -Inf")
      
    }
  }
}

######
# 7.fexp.prime
######
fexp.prime <- function (x, lower) 
{
  if (length(lower) > 1) {
    mapply(fexp.prime, x, lower)
  }
  else {
    if (lower == -Inf) {
      return(1)
    }
    else { 
      stop("\"lower\" should be -Inf")
      
    }
  }
}

######
# 8.tilted.expe
######
tilted.expe <- function (x, tau) 
{
  tabs <- x
  tabs[x >= 0] <- tabs[x >= 0]^2 * tau
  tabs[x < 0] <- tabs[x < 0]^2 * (1- tau)
  tabs
}

######
# 9. tilted.expe.prime
######
tilted.expe.prime <- function (x, tau) 
{
  tabs <- x
  tabs[x >= 0] <- 2*tabs[x >= 0] * tau
  tabs[x < 0] <- 2*tabs[x < 0] * (1- tau)
  tabs
}

######
# 10.ernn.initialize 
######
ernn.initialize <- function (x, y, n.hidden) 
{
  W1 <- matrix(runif((ncol(x) + 1) * n.hidden, -0.5, 0.5), 
               ncol(x) + 1, n.hidden)
  W2 <- matrix(runif((n.hidden + 1) * ncol(y), -0.5, 0.5), 
               n.hidden + 1, ncol(y))
  c(W1, W2)
}

######
# 11.ernn.predict
######
ernn.predict <- function (x, parms) 
{
  if (!is.matrix(x))
    stop("\"x\" must be a matrix")
  weights <- parms$weights
  lower <- parms$lower
  # eps <- min(parms$eps.seq)
  Th <- parms$Th
  x.center <- parms$x.center
  x.scale <- parms$x.scale
  y.center <- parms$y.center
  y.scale <- parms$y.scale
  lower <- (lower - y.center)/y.scale
  x <- sweep(x, 2, x.center, "-")
  x <- sweep(x, 2, x.scale, "/")
  y.bag <- matrix(0, ncol = length(weights), nrow = nrow(x))
  for (i in seq_along(weights)) {
    y.bag[, i] <- ernn.eval(x, weights[[i]]$W1, weights[[i]]$W2, 
                            lower, Th)
    y.bag[, i] <- y.bag[, i] * y.scale + y.center
  }
  y.bag
}

######
# 12.ernn.eval
######
ernn.eval <- function (x, W1, W2, lower, Th)
{
  x <- cbind(x, 1)
  h1 <- x %*% W1
  y1 <- Th(h1)
  aug.y1 <- cbind(y1, 1)
  y2 <- aug.y1 %*% W2
  y2 <- fexp(y2, lower)
  y2
}

######
# 13.huber and huber.prime
######
huber <- function (x, eps) 
{
  h <- ifelse(abs(x) > eps, abs(x) - eps/2, (x^2)/(2 * eps))
  h[is.nan(h)] <- 0
  h
}

huber.prime <- function (x, eps) ## Derivative of the Huber norm function. 
{                               
  dh <- x/eps
  dh[x > eps] <- 1
  dh[x < -eps] <- -1
  dh[is.nan(dh)] <- 0
  dh
}








