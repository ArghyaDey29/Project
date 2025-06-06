propA=NULL
for(j in 1:500)
{
delta=c(rep(1,4),rep(0,4)) #---first 4 patients for treatmentA and Last 4 for treatmentB
CA=NULL
CB=NULL
ta=NULL
tb=NULL
delA=NULL
delB=NULL
RA=NULL
rateA=NULL
rateB=NULL
lamA=NULL
lamB=NULL
phiA=NULL
phiB=NULL
gamA=NULL
gamB=NULL
rhoA=NULL
rhoB=NULL
countA=0
countB=0
for( i in 1:4){
marshal_olkin_sampleA <- function(rateA, lamA, phiA, size=1) {
   x1=rexp(1,rateA)
  x2=rexp(1,lamA)
  x3=rexp(1,phiA)	
  y1 <- min(x1, x3)
  y2 <- min(x2, x3) 
 XA=y1
 CA=y2
  return(cbind(XA, CA))
}
samplesA <- marshal_olkin_sampleA(rateA=9, lamA=6, phiA=0.9, size=1)
	ta[i]=min(samplesA[,1],samplesA[,2])
	tb[i]=0
if(samplesA[,1]<=samplesA[,2])
{
	delA[i]=1;delB[i]=0
}
else
{
		delA[i]=0;delB[i]=0
}
}
for(i in 5:8){
marshal_olkin_sampleB <- function(rateB, lamB, phiB, size=1) {
  x1=rexp(1,rateB)
  x2=rexp(1,lamB)
  x3=rexp(1,phiB)	
  y1 <- min(x1, x3)
  y2 <- min(x2, x3) 

 XB=y1
 CB=y2
  return(cbind(XB, CB))
}
samplesB <- marshal_olkin_sampleB(rateB=4, lamB=5, phiB=0.7, size=1)
	tb[i]=min(samplesB[,1],samplesB[,2])
	ta[i]=0
if(samplesB[,1]<=samplesB[,2])
{
	delB[i]=1;delA[i]=0
}
else
{
		delA[i]=0;delB[i]=0
}
}
countA=4;countB=4
gamA[1]=(sum(delta*delA)+0.5)/(sum(ta)+1)
gamB[1]=(sum((1-delta)*delB)+0.5)/(sum(tb)+1)
lamA[1]=(sum(delta*(1-delA))+0.5)/(sum(ta)+1)
lamB[1]=(sum((1-delta)*(1-delB))+0.5)/(sum(tb)+1)
RA[1]=(gamA[1]+lamA[1]+0.5)/((gamA[1]+gamB[1]+lamA[1]+lamB[1])+1)
RA[1]
gamA[1]
gamB[1]
lamA[1]
lamB[1]
phiA=0.9
phiB=0.7
rhoA[1]=lamA[1]+phiA
rhoB[1]=lamA[1]+phiB
for(i in 9:100)
{
	if(runif(1)<=RA[i-8])
	{
		delta[i]=1
		ta[i]=min(rexp(1,rate=gamA[i-8]),rexp(1,rate=rhoA[i-8]));tb[i]=0
		if(rexp(1,rate=gamA[i-8])<=rexp(1,rate=rhoA[i-8]))
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
		
		tb[i]=min(rexp(1,rate=gamB[i-8]),rexp(1,rate=rhoB[i-8]));ta[i]=0
		if(rexp(1,rate=gamB[i-8])<=rexp(1,rate=rhoB[i-8]))
		{
			delB[i]=1;delA[i]=0

		}
		else
		{
			delA[i]=0;delB[i]=0
		}
		countB=countB+1
	}
gamA[i-7]=(sum(delta*delA)+0.5)/(sum(ta)+1)
gamB[i-7]=(sum((1-delta)*delB)+0.5)/(sum(tb)+1)
lamA[i-7]=(sum(delta*(1-delA))+0.5)/(sum(ta)+1)
lamB[i-7]=(sum((1-delta)*(1-delB))+0.5)/(sum(tb)+1)
RA[i-7]=(gamA[i-7]+lamA[i-7]+0.5)/((gamA[i-7]+gamB[i-7]+lamA[i-7]+lamB[i-7])+1)
rateA[i-7]=gamA[i-7]-phiA
rhoA[i-7]=lamA[i-7]+phiA
rateB[i-7]=gamB[i-7]-phiB
rhoB[i-7]=lamB[i-7]+phiB

	}
propA[j]=countA/100
}
mean_propA=mean(propA)
se_propA=sd(propA)
mean_propA
se_propA
