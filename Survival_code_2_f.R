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
  for(i in 1:8)
  {
  	C_n[i]=rexp(1,rate=4)
  }
  C_n
  for(i in 1:4)
  {    XA[i]=rexp(1,rate=22);XB[i]=0
  	ta[i]=min(XA[i],C_n[i])
  	tb[i]=0
  	if(XA[i]<=C_n[i]) #-- if uncensored then RV takes value 1 
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
  	XB[i]=rexp(1,rate=3);XA[i]=0
  	tb[i]=min(XB[i],C_n[i])
  	ta[i]=0
  	if(XB[i]<=C_n[i])
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
  C_n
  countA=4;countB=4
  #-- estimated value of parameter using MLE

  thetaA[1]=(sum(delta*delA)+0.5)/(sum(ta)+0.5)
  thetaB[1]=(sum((1-delta)*delB)+0.5)/(sum(tb)+0.5)
  lambda[1]=(sum(delta*(1-delA))+sum((1-delta)*(1-delB))+0.5)/(sum(ta)+sum(tb)+0.5)
  RA[1]=(2*thetaA[1]+lambda[1]+0.5)/(2*(thetaA[1]+thetaB[1]+lambda[1])+0.5)
  RA[1]
thetaA[1]
thetaB[1]
lambda[1]
  #-- allocation wrt the Probability obtain by 8 patients and continue update

  for(i in 9:100)
  {	
	  X_A=NULL	
	  X_B=NULL
	  C=NULL
		  X_A[i-8]=rexp(1,rate=thetaA[i-8])
		  X_B[i-8]=rexp(1,rate=thetaB[i-8])
		  C[i-8]=rexp(1,rate=lambda[i-8])
	    if(runif(1)<=RA[i-8])
	    {
		    delta[i]=1
		    ta[i]=min(X_A[i-8],C[i-8]);tb[i]=0
		    XA=c(XA,X_A[i-8]);C_n=c(C_n,C[i-8]);XB=c(XB,0)
		    ta=c(ta,ta[i-8]);tb=c(tb,tb[i])
		    if(X_A[i-8]<=C[i-8])
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
  		tb[i]=min(X_B[i-8],C[i-8]);ta[i]=0
  		XB=c(XB,X_B[i-8]);C_n=c(C_n,C[i-8]);XA=c(XA,0)
  		tb=c(tb,tb[i]);ta=c(ta,ta[i])
  		if(X_B[i-8]<=C[i-8])
  		{
  			delB[i]=1;delA[i]=0
  		}
  		else
  		{
  			delA[i]=0;delB[i]=0
  		}
  		countB=countB+1
  	}
  thetaA[i-7]=(sum(delA)+0.5)/(sum(ta)+0.5)
  thetaB[i-7]=(sum(delB)+0.5)/(sum(tb)+0.5)
  lambda[i-7]=(sum(1-delA)+sum(1-delB)+0.5)/(sum(ta)+sum(tb)+0.5)
  RA[i-7]=(2*thetaA[i-7]+lambda[i-7]+0.5)/(2*(thetaA[i-7]+thetaB[i-7]+lambda[i-7])+0.5)
  }
propA[j]=countA/100
}
mean_propA=mean(propA)
se_propA=sd(propA)
mean_propA
se_propA 
length(ta[-1])
P=as.vector(XA+XB)
Q=cbind(P,C_n)
