### Step 1: Generate empirical critical value under H0 (thetaA = thetaB)
Zn_H0 <- numeric(5000)
for (j in 1:5000) {
  n <- 100
  delta <- c(rep(1, 4), rep(0, 4))
  XA <- XB <- C <- ya <- yb <- delA <- delB <- numeric(n)

  for (i in 1:4) {
    XA[i] <- rexp(1, rate = 6)  # same rate => thetaA = thetaB
    C[i] <- rexp(1, rate = 5)
    ya[i] <- min(XA[i], C[i])
    delA[i] <- as.numeric(XA[i] <= C[i])
  }

  for (i in 5:8) {
    XB[i] <- rexp(1, rate = 6)
    C[i] <- rexp(1, rate = 5)
    yb[i] <- min(XB[i], C[i])
    delB[i] <- as.numeric(XB[i] <= C[i])
  }

  XA_init <- XA[1:4][delA[1:4] == 1]
  XB_init <- XB[5:8][delB[5:8] == 1]
  if (length(XA_init) == 0 || length(XB_init) == 0) next

  thetaA_hat <- length(XA_init) / sum(XA_init)
  thetaB_hat <- length(XB_init) / sum(XB_init)
  theta_hat <- (length(XA_init) + length(XB_init)) / sum(c(XA_init, XB_init))

  pA <- (thetaA_hat + 0.5 * theta_hat) / (thetaA_hat + thetaB_hat + theta_hat)

  for (i in 9:100) {
    C[i] <- rexp(1, rate = 5)
    if (runif(1) <= pA) {
      delta[i] <- 1
      XA[i] <- rexp(1, rate = 10)
      ya[i] <- min(XA[i], C[i])
      delA[i] <- as.numeric(XA[i] <= C[i])
    } else {
      delta[i] <- 0
      XB[i] <- rexp(1, rate = 10)
      yb[i] <- min(XB[i], C[i])
      delB[i] <- as.numeric(XB[i] <= C[i])
    }
   }

  ya_total <- sum(ya)
  yb_total <- sum(yb)
  if (ya_total == 0 || yb_total == 0) next

  theta_An <- sum(delta * delA) / ya_total
  theta_Bn <- sum((1 - delta) * delB) / yb_total
  theta_0_n_H0 <- (sum(delta * (1 - delA)) + sum((1 - delta) * (1 - delB))) / (ya_total + yb_total)
  theta_n_H0 <- sum((delta * (1 - delA) * (1 - delB)) + ((1 - delta) * (1 - delA) * (1 - delB))) / (ya_total + yb_total)

  Zn <- sqrt(n) * (theta_An - theta_Bn) / (2 * sqrt(theta_0_n_H0 * (theta_n_H0 + theta_0_n_H0)))
  Zn_H0[j] <- Zn
}

# Empirical critical value (95th percentile)
crit_val_empirical <- quantile(Zn_H0, probs = 0.95, na.rm = TRUE)
cat("Empirical 95% critical value under H0:", crit_val_empirical, "\n")

### Step 2: Power calculation using empirical cutoff
set.seed(123)
the_val=6:15
po=NULL
Zn_values <- numeric(5000)
rejection <- numeric(5000)
for(k in 1:length(the_val)){
for (j in 1:5000) {
  n <- 100
  delta <- c(rep(1, 4), rep(0, 4))
  XA <- XB <- C <- ya <- yb <- delA <- delB <- numeric(n)

  for (i in 1:4) {
    XA[i] <- rexp(1, rate = the_val[k])  # thetaA under H1
    C[i]  <- rexp(1, rate = 5)
    ya[i] <- min(XA[i], C[i])
    delA[i] <- as.numeric(XA[i] <= C[i])
  }

  for (i in 5:8) {
    XB[i] <- rexp(1, rate = 6)   # thetaB
    C[i]  <- rexp(1, rate = 5)
    yb[i] <- min(XB[i], C[i])
    delB[i] <- as.numeric(XB[i] <= C[i])
  }

  XA_init <- XA[1:4][delA[1:4] == 1]
  XB_init <- XB[5:8][delB[5:8] == 1]
  if (length(XA_init) == 0 || length(XB_init) == 0) next

  thetaA_hat <- length(XA_init) / sum(XA_init)
  thetaB_hat <- length(XB_init) / sum(XB_init)
  theta_hat <- (length(XA_init) + length(XB_init)) / sum(c(XA_init, XB_init))

  pA <- (thetaA_hat + 0.5 * theta_hat) / (thetaA_hat + thetaB_hat + theta_hat)

  for (i in 9:100) {
    C[i] <- rexp(1, rate = 5)
    if (runif(1) <= pA) {
      delta[i] <- 1
      XA[i] <- rexp(1, rate = the_val[k])
      ya[i] <- min(XA[i], C[i])
      delA[i] <- as.numeric(XA[i] <= C[i])
    } else {
      delta[i] <- 0
      XB[i] <- rexp(1, rate = 6)
      yb[i] <- min(XB[i], C[i])
      delB[i] <- as.numeric(XB[i] <= C[i])
    }
  }

  ya_total <- sum(ya)
  yb_total <- sum(yb)
  if (ya_total == 0 || yb_total == 0) next

  theta_An <- sum(delta * delA) / ya_total
  theta_Bn <- sum((1 - delta) * delB) / yb_total
  theta_0_n_H0 <- (sum(delta * (1 - delA)) + sum((1 - delta) * (1 - delB))) / (ya_total + yb_total)
  theta_n_H0 <- sum((delta * (1 - delA) * (1 - delB)) + ((1 - delta) * (1 - delA) * (1 - delB))) / (ya_total + yb_total)

  Zn <- sqrt(n) * (theta_An - theta_Bn) / (2 * sqrt(theta_0_n_H0 * (theta_n_H0 + theta_0_n_H0)))
  Zn_values[j] <- Zn
  rejection[j] <- as.numeric(Zn > crit_val_empirical)
}
po[k]=sum(rejection, na.rm = TRUE) / 5000
}

cat("ThetaA values:\n")
print(the_val)
cat("Corresponding powers by empirical cut off (adaptive):\n")
print(po)

##
set.seed(123)
the_val=6:15
po2=NULL
Zn_values <- numeric(5000)
rejection <- numeric(5000)
for(k in 1:length(the_val)){
for (j in 1:5000) {
  n <- 100
  delta <- c(rep(1, 4), rep(0, 4))
  XA <- XB <- C <- ya <- yb <- delA <- delB <- numeric(n)

  for (i in 1:4) {
    XA[i] <- rexp(1, rate = the_val[k])  # thetaA under H1
    C[i]  <- rexp(1, rate = 5)
    ya[i] <- min(XA[i], C[i])
    delA[i] <- as.numeric(XA[i] <= C[i])
  }

  for (i in 5:8) {
    XB[i] <- rexp(1, rate = 6)   # thetaB
    C[i]  <- rexp(1, rate = 5)
    yb[i] <- min(XB[i], C[i])
    delB[i] <- as.numeric(XB[i] <= C[i])
  }

  XA_init <- XA[1:4][delA[1:4] == 1]
  XB_init <- XB[5:8][delB[5:8] == 1]
  if (length(XA_init) == 0 || length(XB_init) == 0) next

  thetaA_hat <- length(XA_init) / sum(XA_init)
  thetaB_hat <- length(XB_init) / sum(XB_init)
  theta_hat <- (length(XA_init) + length(XB_init)) / sum(c(XA_init, XB_init))

  pA <- (thetaA_hat + 0.5 * theta_hat) / (thetaA_hat + thetaB_hat + theta_hat)

  for (i in 9:100) {
    C[i] <- rexp(1, rate = 5)
    if (runif(1) <= pA) {
      delta[i] <- 1
      XA[i] <- rexp(1, rate = the_val[k])
      ya[i] <- min(XA[i], C[i])
      delA[i] <- as.numeric(XA[i] <= C[i])
    } else {
      delta[i] <- 0
      XB[i] <- rexp(1, rate = 6)
      yb[i] <- min(XB[i], C[i])
      delB[i] <- as.numeric(XB[i] <= C[i])
    }
  }

  ya_total <- sum(ya)
  yb_total <- sum(yb)
  if (ya_total == 0 || yb_total == 0) next

  theta_An <- sum(delta * delA) / ya_total
  theta_Bn <- sum((1 - delta) * delB) / yb_total
  theta_0_n_H0 <- (sum(delta * (1 - delA)) + sum((1 - delta) * (1 - delB))) / (ya_total + yb_total)
  theta_n_H0 <- sum((delta * (1 - delA) * (1 - delB)) + ((1 - delta) * (1 - delA) * (1 - delB))) / (ya_total + yb_total)

  Zn <- sqrt(n) * (theta_An - theta_Bn) / (2 * sqrt(theta_0_n_H0 * (theta_n_H0 + theta_0_n_H0)))
  Zn_values[j] <- Zn
  rejection[j] <- as.numeric(Zn > 1.645)
}
po2[k]=sum(rejection, na.rm = TRUE) / 5000
}
cat("ThetaA values:\n")
print(the_val)
cat("Corresponding powers by normal dist cut off (adaptive):\n")
print(po2)
par(mfrow=c(1,2))
plot(the_val,po,main="Power Curve\n (By Empirical Cut off)",type="l",lwd=3,xlab="Values of Theta",ylab="Power")
plot(the_val,po2,main="Power Curve\n (By Normal Dist Cut off)",type="l",lwd=3,xlab="Values of Theta",ylab="Power")