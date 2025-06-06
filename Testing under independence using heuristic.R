set.seed(123)

generate_Tn_adaptive <- function(n, thetaA, thetaB, theta) {
  n_initial <- 8
  delta <- numeric(n)     # treatment indicator: 1 for A, 0 for B
  delta[1:4] <- 1
  delta[5:8] <- 0

  XA <- rexp(n, rate = thetaA)
  XB <- rexp(n, rate = thetaB)
  C  <- rexp(n, rate = theta)

  time <- numeric(n)
  deltaA <- numeric(n)    # uncensoring indicator for treatment A
  deltaB <- numeric(n)    # uncensoring indicator for treatment B

  # Initial fixed allocation for first 8 patients
  for (i in 1:n_initial) {
    if (delta[i] == 1) {
      time[i] <- min(XA[i], C[i])
      deltaA[i] <- as.numeric(XA[i] <= C[i])
    } else {
      time[i] <- min(XB[i], C[i])
      deltaB[i] <- as.numeric(XB[i] <= C[i])
    }
  }

  # Estimate hazard rates from initial uncensored data
  XA_init <- XA[1:4][deltaA[1:4] == 1]
  XB_init <- XB[5:8][deltaB[5:8] == 1]

  if (length(XA_init) == 0 || length(XB_init) == 0) return(NA)

  thetaA_hat <- length(XA_init) / sum(XA_init)
  thetaB_hat <- length(XB_init) / sum(XB_init)
  theta_hat  <- (length(XA_init) + length(XB_init)) / sum(c(XA_init, XB_init))

  # Updated adaptive allocation
  pA <- (thetaA_hat + 0.5 * theta_hat) / (thetaA_hat + thetaB_hat + theta_hat)

  for (i in (n_initial + 1):n) {
    delta[i] <- rbinom(1, 1, pA)
  }

  # Simulate remaining patients
  for (i in (n_initial + 1):n) {
    if (delta[i] == 1) {
      time[i] <- min(XA[i], C[i])
      deltaA[i] <- as.numeric(XA[i] <= C[i])
    } else {
      time[i] <- min(XB[i], C[i])
      deltaB[i] <- as.numeric(XB[i] <= C[i])
    }
  }

  # Compute test statistic
  N_An <- sum(delta * deltaA)
  N_Bn <- sum((1 - delta) * deltaB)

  sum_XA <- sum(delta * XA)
  sum_XB <- sum((1 - delta) * XB)

  if (N_An == 0 || N_Bn == 0 || sum_XA == 0 || sum_XB == 0) return(NA)

  theta_An_hat <- N_An / sum_XA
  theta_Bn_hat <- N_Bn / sum_XB
  theta0n_hat  <- n / ((N_An / theta_An_hat) + (N_Bn / theta_Bn_hat))

  Tn <- sqrt(n) * (theta_An_hat - theta_Bn_hat) / (sqrt(8) * theta0n_hat)

  return(Tn)
}

# ==== Simulation Settings ====
n <- 100
thetaB <- 6
theta0 <- 6
theta <- 4
n_iter <- 10000

# ==== Simulate under H0 ====
Tn_H0 <- numeric(n_iter)
for (j in 1:n_iter) {
  Tn_H0[j] <- generate_Tn_adaptive(n, theta0, theta0, theta)
}
Tn_H0 <- Tn_H0[!is.na(Tn_H0)]
critical_val <- quantile(Tn_H0, 0.95)
cat("Critical value:", critical_val, "\n")

# ==== Simulate under H1 (thetaA = 12) ====
Tn_H1 <- numeric(n_iter)
for (j in 1:n_iter) {
  Tn_H1[j] <- generate_Tn_adaptive(n, 12, thetaB, theta)
}
Tn_H1 <- Tn_H1[!is.na(Tn_H1)]

# ==== Empirical Power Calculation ====
n_exceeds <- sum(Tn_H1 > critical_val, na.rm = TRUE)
emp_power <- n_exceeds / n_iter
cat("Empirical Power at thetaA = 12 (adaptive):", emp_power, "\n")


# ==== Power Curve for thetaA from 6 to 15 ====
thetaA_vals <- 6:15
power_vals <- numeric(length(thetaA_vals))

for (i in seq_along(thetaA_vals)) {
  temp_Tn <- numeric(n_iter)
  for (j in 1:n_iter) {
    temp_Tn[j] <- generate_Tn_adaptive(n, thetaA_vals[i], thetaB, theta)
  }
  temp_Tn <- temp_Tn[!is.na(temp_Tn)]
  power_vals[i] <- mean(temp_Tn > critical_val)
}

cat("ThetaA values:\n")
print(thetaA_vals)
cat("Corresponding powers (adaptive):\n")
print(power_vals)
##
thetaA_vals <- 6:15
power_vals2 <- numeric(length(thetaA_vals))

for (i in seq_along(thetaA_vals)) {
  temp_Tn <- numeric(n_iter)
  for (j in 1:n_iter) {
    temp_Tn[j] <- generate_Tn_adaptive(n, thetaA_vals[i], thetaB, theta)
  }
  temp_Tn <- temp_Tn[!is.na(temp_Tn)]
  power_vals2[i] <- mean(temp_Tn > 1.625)
}

cat("ThetaA values:\n")
print(thetaA_vals)
cat("Corresponding powers (adaptive):\n")
print(power_vals2)
# ==== Plotting ====
par(mfrow = c(1,2))

plot(thetaA_vals, power_vals, type = "l", pch = 19,lwd = 2,main = "Power Curve (Adaptive)", xlab = "ThetaA", ylab = "Empirical Power")
plot(thetaA_vals, power_vals2, type = "l", pch = 19, lwd = 2,main = "Power Curve(Adaptive)\n with Normal dist Cutoff", xlab = "ThetaA", ylab = "Empirical Power")