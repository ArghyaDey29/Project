propA=NULL
for(j in 1:5000)
{
  delta=c(rep(1,4),rep(0,4)) #---first 4 patients for tr_A and Last 4 for tr_B
  C_n=NULL
  XA=NULL
  XB=NULL
  ta=NULL
  tb=NULL
  delA=NULL #-- indicator RV if the patient allocated to tr_A Censored or not
  delB=NULL 
  RA=NULL
  thetaA=NULL
  thetaB=NULL
  lambda=NULL
  countA=0
  countB=0
  nu=4
  C_A=NULL
  C_B=NULL
  for(i in 1:4)
  {    XA[i]=rexp(1,rate=20);XB[i]=0
    C_A[i]=rexp(1,rate=nu*7)
    ta[i]=min(XA[i],C_A[i])
  tb[i]=0
  if(XA[i]<=C_A[i]) #-- if uncensored then RV takes value 1 
  {
    delA[i]=1;delB[i]=0 
  }
  else
  {
    delA[i]=0;delB[i]=0
  }
  }
  for(i in 5:8)
  {
    XB[i]=rexp(1,rate=4);XA[i]=0
    C_B[i]=rexp(1,rate=nu*4)
    tb[i]=min(XB[i],C_B[i])
    ta[i]=0
    if(XB[i]<=C_B[i])
    {
      delA[i]=0;delB[i]=1
    }
    else
    {
      delA[i]=0;delB[i]=0
    }
  }
  ta
  tb
  XA
  XB
  countA=4;countB=4
  #-- estimated value of parameter using MLE
  
  thetaA[1]=((sum(delta)+0.5)/((sum(ta)*(1+nu))+0.5))
  thetaB[1]=((sum((1-delta))+0.5)/((sum(tb)*(1+nu))+0.5))
  RA[1]=thetaA[1]/(thetaA[1]+thetaB[1])
  RA[1]
  thetaA[1]
  thetaB[1]
  #-- allocation wrt the Probability obtain by 8 patients and continue update
  
  for(i in 9:100)
  {	
    X_A=NULL	
    X_B=NULL
    CA=NULL
    CB=NULL
    X_A[i-8]=rexp(1,rate=thetaA[i-8])
    X_B[i-8]=rexp(1,rate=thetaB[i-8])
    CA[i-8]=rexp(1,rate=nu*thetaA[i-8])
    CB[i-8]=rexp(1,rate=nu*thetaB[i-8])
    if(runif(1)<=RA[i-8])
    {
      delta[i]=1
      ta[i]=min(X_A[i-8],CA[i-8]);tb[i]=0
      XA=c(XA,X_A[i-8]);C_A=c(C_A,CA[i-8]);XB=c(XB,0)
      ta=c(ta,ta[i-8]);tb=c(tb,tb[i])
      if(X_A[i-8]<=CA[i-8])
      {
        delA[i]=1;delB[i]=0
      }
      else
      {
        delA[i]=0;delB[i]=0
      }
      countA=countA+1
    }
    else
    {
      delta[i]=0
      tb[i]=min(X_B[i-8],CB[i-8]);ta[i]=0
      XB=c(XB,X_B[i-8]);C_B=c(C_B,CB[i-8]);XA=c(XA,0)
      tb=c(tb,tb[i]);ta=c(ta,ta[i])
      if(X_B[i-8]<=CB[i-8])
      {
        delB[i]=1;delA[i]=0
      }
      else
      {
        delA[i]=0;delB[i]=0
      }
      countB=countB+1
    }
    thetaA[i-7]=(sum(delA)+0.5)/((sum(ta)*(1+nu))+0.5)
    thetaB[i-7]=(sum(delB)+0.5)/((sum(tb)*(1+nu))+0.5)
    RA[i-7]=thetaA[i-7]/(thetaA[i-7]+thetaB[i-7])
  }
  propA[j]=countA/100
}
mean_propA=mean(propA)
se_propA=sd(propA)
mean_propA
se_propA 
length(ta[-1])
P=as.vector(XA+XB)
Q=as.vector(CA+CB)
