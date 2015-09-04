################################################################################################################
##
## Price simulations
##
################################################################################################################
#################simulate GARCH: the error of model 
library(quadprog)
set.seed(1288)
omega=rnorm(100,mean=0,sd=5)^2
alpha=runif(100,min=0.1,max=0.3)
beta=runif(100,min=0.5,max=0.7)


## A simple GARCH(1,1) simulator
garch11=function(omega, alpha, beta){
  sig=c()
  x=c()
  sig[1]=0
  x[1]=0
  for(i in 2:400){
    sig[i]=sqrt(omega+alpha*x[i-1]^2+beta*sig[i-1]^2)
    x[i]=sig[i]*rnorm(1,0,1)}
  return(cbind(sig[151:400],x[151:400])) #burn in 
}

sig=matrix(0,nrow=250,ncol=100) # the sd of error
z=matrix(0,nrow=250,ncol=100) # error
for(i in 1:100){
  ans = garch11(omega=omega[i],alpha = alpha[i], beta = beta[i])
  sig[,i]=ans[,1]
  z[,i]=ans[,2]
}



##############simulate the coefficient A
A=matrix(0,nrow=250,ncol=100)
for(i in 1:250){
  A[i,]=rnorm(100,1,0.0005)
}

############simulate the price!
price=matrix(0,nrow=251,ncol=100)
price[1,]=runif(100,min=100,max=300) #initial price

for(i in 1:250){
  price[i+1,]=price[i,]*A[i,]+z[i,]}

########burn in
price=price[101:251,]
A=A[101:250,]
sig=sig[101:250,]
##########

for(i in 1:10){
  plot(price[,i],type='l')
}
################################################################################################################
##
## initial optimizer
##
################################################################################################################

opt1<-function(R,C,price,lb,ub,gamma){
  l=length(R);
  Amat <- cbind(matrix(price,nrow=l,ncol=1),diag(price),-diag(price));
  bvec<-c(10000,rep(lb*10000,l),rep(-ub*10000,l));
  k<-solve.QP(C,0.5/gamma*R,Amat,bvec,meq=1)
  return(k)
}




################################################################################################################
##
##  optimizer with transaction cost constraint 
##
##
##  E=expectation of terminal stock price
##  CV= covariance variance matrix of stock price
##  IH= initial houlding of each stocks
##  lb,ub= lower and uper bond of weight eg.-0.05~0.05
##  TC= transaction cost eg. 1% of value change
##
################################################################################################################



OPT=function(E,CV,price,IH,lb,ub,TC,gamma){
  l=length(E)
  
  Dmat=diag(c(diag(CV),rep(0.1^10,100*2)))
  dvec=matrix(c(E,rep(0,100*2)),ncol=100*3,nrow=1)
  
  tc=TC # transaction cost: such as 1% of value change
  
  #constraint 1, capital 
  A1=c(price,price*tc,price*(-tc))
  B1=IH%*%price
  
  #constraint 2, # of initial shares= #shares-#buy+#sell
  A2=rbind(diag(l),-diag(l),-diag(l))
  B2=IH
  
  #constraint 3, Stocks buy  >0
  A3=rbind(matrix(0,nrow=l,ncol=l),diag(l),matrix(0,nrow=l,ncol=l))
  B3=rep(0,l)
  
  #constraint 4, Stocks sell <0
  A4=rbind(matrix(0,nrow=l,ncol=l),matrix(0,nrow=l,ncol=l),-diag(l))
  B4=rep(0,l)
  
  #constraint 5, initial capital>= adjusted capital
  A5=-1*price
  B5=-1*IH%*%price
  
  #constraint 6 lower bound, weight bound 
  A6=rbind(diag(price),matrix(0,nrow=l,ncol=l),matrix(0,nrow=l,ncol=l))
  B6=rep(lb*price%*%IH,l)
  
  #constraint 7 uper bound, weight bound
  A7=rbind(diag(price),matrix(0,nrow=l,ncol=l),matrix(0,nrow=l,ncol=l))
  B7=rep(-ub*price%*%IH,l)
  
  
  Amat=cbind(A1,A2,A3,A4,A5,A6,A7)
  bvec=c(B1,B2,B3,B4,B5,B6,B7)
  
  ans=solve.QP(Dmat, 0.5/gamma*dvec, Amat, bvec, meq=l+1)$solution
  
  solution=list(portfolio = ans[1:l], buy=ans[(l+1):(2*l)], sell=ans[(2*l+1):(3*l)])
  return(solution)
}


################################################################################################################
##
## find best porfolios from T0-T150
##  GAMMA=1
## maxDrawdown(R)# 0.001656589
## SharpeRatio(R,Rf) # 0.3246886
## InformationRatio(R,Rf)# 1.124881
## AverageDrawdown(R) #0.0005390073
##
################################################################################################################
library(quadprog)
C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
turnover=matrix(0,nrow=150,ncol=1)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=1)$solution
turnover[1,]=0
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  ans<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.01,gamma=1)
  k[i,]<-ans$portfolio
  turnover[i,]<-(ans$buy-ans$sell)%*%price[i,]/(k[i,]%*%price[i,])
}



##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


maxDrawdown(R)# 0.001656589
SharpeRatio(R,Rf) # 0.3246886
InformationRatio(R,Rf)# 1.124881
AverageDrawdown(R) #0.0005390073

################################################################################################################
##
## find best porfolios from T0-T149
##  GAMMA=0.1
##  maxDrawdown(R)#0.001710163
##  SharpeRatio(R,Rf,FUN = "StdDev")#0.7679199
##  InformationRatio(R,Rf)#2.666087
##  AverageDrawdown(R) #0.0003653663
################################################################################################################



C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=0.1)$solution
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  k[i,]<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.01,gamma=0.1)$portfolio
}


##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


maxDrawdown(R)#0.001710163
SharpeRatio(R,Rf,FUN = "StdDev")#0.7679199
InformationRatio(R,Rf)#2.666087
AverageDrawdown(R) #0.0003653663

################################################################################################################
##
## find best porfolios from T0-T149
##  GAMMA=0.01
## maxDrawdown(R)#0.007756667
## SharpeRatio(R,Rf,FUN = "StdDev")#0.703222
## InformationRatio(R,Rf)#2.44937
## AverageDrawdown(R) #0.001434193
##
################################################################################################################



C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=0.01)$solution
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  k[i,]<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.01,gamma=0.01)$portfolio
}


##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


maxDrawdown(R)#0.007756667
SharpeRatio(R,Rf,FUN = "StdDev")#0.703222
InformationRatio(R,Rf)#2.44937
AverageDrawdown(R) #0.001434193


################################################################################################################
##
## find best porfolios from T0-T149
##  GAMMA=0.05
##  maxDrawdown(R)#0.001835507
## SharpeRatio(R,Rf,FUN = "StdDev")#0.7917611
## InformationRatio(R,Rf)#2.750856
## AverageDrawdown(R) #0.0005482151
################################################################################################################



C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=0.05)$solution
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  k[i,]<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.01,gamma=0.05)$portfolio
}


##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


 maxDrawdown(R)#0.001835507
 SharpeRatio(R,Rf,FUN = "StdDev")#0.7917611
 InformationRatio(R,Rf)#2.750856
 AverageDrawdown(R) #0.0005482151







################################################################################################################
##
## find best porfolios from T0-T149
##  GAMMA=1*10000/X   Better than gamma=1 !
##  maxDrawdown(R)#0.001653762
##  SharpeRatio(R,Rf,FUN = "StdDev")#0.326596
##  InformationRatio(R,Rf)#1.131499
##  AverageDrawdown(R) #0.0005367311
##
################################################################################################################



C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=1)$solution
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  gamma=10000/(as.numeric(k[i-1,]%*%price[i,]))
  k[i,]<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.01,gamma)$portfolio
}


##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


maxDrawdown(R)#0.001653762
SharpeRatio(R,Rf,FUN = "StdDev")#0.326596
InformationRatio(R,Rf)#1.131499
AverageDrawdown(R) #0.0005367311


################################################################################################################
##
## find best porfolios from T0-T149
##  GAMMA=0.1*10000/X   Better than gamma=0.1 !
##  maxDrawdown(R)#0.001709914
##  SharpeRatio(R,Rf,FUN = "StdDev")#0.7707199
##  InformationRatio(R,Rf)#2.675914
##  AverageDrawdown(R) #0.0003630485
##
################################################################################################################



C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=0.1)$solution
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  gamma=0.1*10000/(as.numeric(k[i-1,]%*%price[i,]))
  k[i,]<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.01,gamma)$portfolio
}


##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


maxDrawdown(R)#0.001709914
SharpeRatio(R,Rf,FUN = "StdDev")#0.7707199
InformationRatio(R,Rf)#2.675914
AverageDrawdown(R) #0.0003630485




################################################################################################################
##
## find best porfolios from T0-T149
##  TC=0.02 P&L slightly worse than 0.01
##  GAMMA=0.1
##  maxDrawdown(R)#0.001710163
##  SharpeRatio(R,Rf,FUN = "StdDev")#0.7679199
##  InformationRatio(R,Rf)#2.666087
##  AverageDrawdown(R) #0.0003653663
################################################################################################################



C1=diag(sig[1,]^2)
mu1=matrix(c(A[1,]*price[1,]),ncol=100,nrow=1)

k=matrix(0,nrow=150,ncol=100)
k[1,]=opt1(mu1,C1,price[1,],lb=-0.05,ub=0.05,gamma=0.1)$solution
for(i in 2:150){
  C=diag(sig[i,]^2)
  mu=matrix(c(A[i,]*price[i,]),ncol=100,nrow=1)
  k[i,]<-OPT(mu,C,price[i,],k[i-1,],-0.05,0.05,TC=0.02,gamma=0.1)$portfolio
}


##pay and loss
pl=c()
for(i in 1:150){
  pl[i]=price[i+1,]%*%k[i,]-10000
}
plot(c(0,pl[1:150]),type='l')
abline(h=0,lty=2)



X=c()
for(i in 1:150){
  X[i]=price[i+1,]%*%k[i,]
}
X=c(10000,X)
Xreturn=diff(X)/X[1:150]
library(PerformanceAnalytics)

x <- as.yearmon(2000 + seq(0, 149)/12)
R=data.frame(Xreturn,row.names=as.Date(x))
Rf=data.frame(rep(0,150),row.names=as.Date(x))
plot(Xreturn,type='l')


maxDrawdown(R)#0.00170964
SharpeRatio(R,Rf,FUN = "StdDev")#0.767639
InformationRatio(R,Rf)#2.665109
AverageDrawdown(R) #0.0003655768


