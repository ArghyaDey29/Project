propA=NULL
n_super_iter <- 5000
for (j in 1:n_super_iter) 
{
  ProbA=NULL
  thetaA=NULL
  thetaB=NULL
  lam=NULL
  XA=NULL
  XB=NULL
  CA=NULL
  X=NULL
  Sam=NULL
  CountA <- 0
  CountB <- 0
  n_inner_iter <- 100 
  #-- equal allocation
  
  SamA <- rexp(n = 4, rate =20) #-- collecting 4 exact response tr_A
  SamB <- rexp(n = 4, rate = 3) #-- collecting 4 exact response for tr_B
  lamA <- rexp(n = 8, rate = 4) #-- 8 censor value of the 8 individual
  
  #-- find the estimated value of parameter using method of moments
  thetaA[1] <- 1 / mean(SamA) 
  thetaB[1] <- 1 / mean(SamB)
  lam[1] <- 1/mean(lamA)
  
  #-- find the probability that tr_A be a better treatment than tr_B 
  ProbA[1] <- (2 * thetaA[1] + lam[1]+0.5) / ((thetaA[1] + thetaB[1] + lam[1])*2+1)
  ProbA[1] #-- based on first 8 samples
  
  for (i in 9:n_inner_iter) 
  {
    runif_sample <- runif(1)    
    if (runif_sample <= ProbA[i-8]) 
    {
      # Case for A (allocated the patient tr_A)
      XA[i-8]=  rexp(1, thetaA[i-8]) #-- exact value
      CA[i-8]=  rexp(1, lam[i-8]) #-- censoring value
      X[i-8] <- min(XA[i-8],CA[i-8]) #-- response
      #-- updated values
      
      SamA=c(SamA,X[i-8]) 
        lamA=c(lamA,X[i-8])
      CountA <- CountA + 1 
    } else {
      # Case for B
      XB[i-8]=    rexp(1, thetaB[i-8])
      CA[i-8]=   rexp(1, lam[i-8])
      X[i-8] <- min(XB[i-8],CA[i-8])
        SamB=c(SamB,X[i-8])
        lamA=c(lamA,X[i-8])
      CountB <- CountB + 1
    }
    #-- updated Parameter and Probability after every allocation
    
    thetaA[i-7] <- 1 / mean(SamA)
    thetaB[i-7] <- 1 / mean(SamB)
    lam[i-7] <- 1/mean(lamA)
    ProbA[i-7] <- (2 * thetaA[i-7] + lam[i-7]+0.5) / ((thetaA[i-7] + thetaB[i-7] + lam[i-7])*2+1)
    
  }
  propA[j]=CountA/100
}
table(propA)
mean_propA=mean(propA)
se_propA=sd(propA)
mean_propA
se_propA
