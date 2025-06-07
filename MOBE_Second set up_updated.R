# ----- Marshall-Olkin Sample Generator -----
marshal_olkin_sample <- function(rate=1, lam=1, phi=1) {
  x1 <- rexp(1, rate)
  x2 <- rexp(1, lam)
  x3 <- rexp(1, phi)	
  y1 <- min(x1, x3)
  y2 <- min(x2, x3) 
  return(c(y1, y2))
}

# ----- Function to Generate d1, d2, correlation -----
generate_D <- function(thetaA=6, thetaB=4, lambA=4, lambB=5, phi_A=0.7, phi_B=0.9) {
  delta <- c(rep(1, 4), rep(0, 4))
  ta <- tb <- delA <- delB <- numeric(100)
  gamA <- gamB <- lamA <- lamB <- RA <- rateA <- rateB <- rhoA <- rhoB <- theta <- lam <- d1 <- d2 <- numeric(100)
  phi <- 0.7

  for (i in 1:4) {
    samplesA <- marshal_olkin_sample(rate=thetaA, lam=lambA, phi=phi_A)
    ta[i] <- min(samplesA[1], samplesA[2])
    delA[i] <- as.integer(samplesA[1] <= samplesA[2])
  }

  for (i in 5:8) {
    samplesB <- marshal_olkin_sample(rate=thetaB, lam=lambB, phi=phi_B)
    tb[i] <- min(samplesB[1], samplesB[2])
    delB[i] <- as.integer(samplesB[1] <= samplesB[2])
  }

  gamA[1] <- (sum(delta * delA) + 0.5) / (sum(ta) + 1)
  gamB[1] <- (sum((1 - delta) * delB) + 0.5) / (sum(tb) + 1)
  lamA[1] <- (sum(delta * (1 - delA)) + 0.5) / (sum(ta) + 1)
  lamB[1] <- (sum((1 - delta) * (1 - delB)) + 0.5) / (sum(tb) + 1)
  RA[1] <- (gamA[1] + lamA[1] + 0.5) / (gamA[1] + gamB[1] + lamA[1] + lamB[1] + 1)
  rateA[1] <- gamA[1] - phi_A
  rateB[1] <- gamB[1] - phi_B
  theta[1] <- (sum(delA * delta) + sum((1 - delta) * delB)) / (sum(ta) + sum(tb)) - phi
  rhoA[1] <- lamA[1] + phi_A
  rhoB[1] <- lamB[1] + phi_B
  lam[1] <- (sum(delta * (1 - delA)) + sum((1 - delta) * (1 - delB))) / (sum(ta) + sum(tb))
  d1[1] <- (rateA[1] - rateB[1]) * sqrt(8) / (2 * (theta[1] + phi))
  d2[1] <- (lamA[1] - lamB[1]) * sqrt(8) / (2 * (theta[1] + phi))

  for (i in 9:100) {
    if (runif(1) <= RA[i - 8]) {
      delta[i] <- 1
      t1 <- rexp(1, gamA[i - 8])
      t2 <- rexp(1, rhoA[i - 8])
      ta[i] <- min(t1, t2); tb[i] <- 0
      delA[i] <- as.integer(t1 <= t2)
    } else {
      delta[i] <- 0
      t1 <- rexp(1, gamB[i - 8])
      t2 <- rexp(1, rhoB[i - 8])
      tb[i] <- min(t1, t2); ta[i] <- 0
      delB[i] <- as.integer(t1 <= t2)
    }

    gamA[i - 7] <- (sum(delta * delA) + 0.5) / (sum(ta) + 1)
    gamB[i - 7] <- (sum((1 - delta) * delB) + 0.5) / (sum(tb) + 1)
    lamA[i - 7] <- (sum(delta * (1 - delA)) + 0.5) / (sum(ta) + 1)
    lamB[i - 7] <- (sum((1 - delta) * (1 - delB)) + 0.5) / (sum(tb) + 1)
    RA[i - 7] <- (gamA[i - 7] + lamA[i - 7] + 0.5) / (gamA[i - 7] + gamB[i - 7] + lamA[i - 7] + lamB[i - 7] + 1)
    rateA[i - 7] <- gamA[i - 7] - phi_A
    rateB[i - 7] <- gamB[i - 7] - phi_B
    rhoA[i - 7] <- lamA[i - 7] + phi_A
    rhoB[i - 7] <- lamB[i - 7] + phi_B
    theta[i - 7] <- (sum(delA * delta) + sum((1 - delta) * delB)) / (sum(ta) + sum(tb)) - phi
    lam[i - 7] <- (sum(delta * (1 - delA)) + sum((1 - delta) * (1 - delB))) / (sum(ta) + sum(tb))
    d1[i - 7] <- (rateA[i - 7] - rateB[i - 7]) * sqrt(i + 8) / (2 * (theta[i - 7] + phi))
    d2[i - 7] <- (lamA[i - 7] - lamB[i - 7]) * sqrt(i + 8) / (2 * (theta[i - 7] + phi))
  }

  r <- cor(d1, d2)
  r <- max(min(r, 0.99), -0.99)  # avoid instability
  return(c(d1[93], d2[93], r))
}

# ----- Function to Calculate Q Statistic -----
generate_Q <- function(thetaA=6, thetaB=4, lambA=4, lambB=5, phi_A=0.7, phi_B=0.9) {
  d <- generate_D(thetaA, thetaB, lambA, lambB, phi_A, phi_B)
  q1 <- sqrt(1 / (1 - d[3]^2)) * sqrt(d[1]^2 + d[2]^2 - 2 * d[3] * d[1] * d[2])
  q2 <- sqrt(1 / (1 - d[3]^2)) * (d[1] - d[3] * d[2])
  q3 <- sqrt(1 / (1 - d[3]^2)) * (d[2] - d[3] * d[1])
  if (d[1] > 0 & d[2] > 0) return(q1)
  else if (d[1] >= d[2] & d[2] < 0) return(q2)
  else if (d[2] >= d[1] & d[1] < 0) return(q3)
  else return(0)
}

# ----- Power Estimation Using Fixed Critical Value -----
set.seed(123)
critical_value <- 1.645  # one-sided 5% level

Q_H1 <- replicate(1000, generate_Q(thetaA = 10, thetaB = 4, lambA = 4, lambB = 4, phi_A = 0.6, phi_B = 0.6))
empirical_power <- mean(Q_H1 > critical_value)

cat("Critical Value Used:", critical_value, "\n")
cat("Empirical Power:", empirical_power, "\n")
set.seed(123)
val=6:18
pow=NULL
options(warn=-1)
for(i in 1:length(val)){
Q_H1 <- replicate(1000, generate_Q(thetaA = val[i], thetaB = 4, lambA = 4, lambB = 4, phi_A = 0.6, phi_B = 0.6))
pow[i]<- mean(Q_H1 > critical_value)
}
pow