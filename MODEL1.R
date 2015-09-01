PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
DATA.T1<-PRICE100[2817:3116,11:30] # Data t-1
DATA.T<-PRICE100[2818:3117,11:30] # Data t
price=as.numeric(DATA.T1[1,])


#######################################################################################################################
##
## fit model X_t=AX_t-1+tau N(0,1)
##
#######################################################################################################################
A=matrix(0,nrow=1,ncol=20)
sd=matrix(0,nrow=1,ncol=20)


  for(i in 1:20){
    y=DATA.T[,i]
    x=DATA.T1[,i]
    L.Model=lm(y~x-1)
    A[1,i]=as.numeric(L.Model$coefficients)
    sd[1,i]=sd(L.Model$residuals)}

#######################################################################################################################
##
## simulations
##
#######################################################################################################################
i=1



s=matrix(0,nrow=301,ncol=5000*20) # 3min run time
for(i in 1:20){
    for(m in 1:5000){
      s[1,(i-1)*5000+m]=price[i] # initial price
      for(n in 1:300){
        s[n+1,(i-1)*5000+m]=A[i]*s[n,(i-1)*5000+m]+rnorm(1,0,sd[i])
      }
    }}

E=matrix(0,nrow=1,ncol=20)
VA=matrix(0,nrow=1,ncol=20)
K=matrix(0,nrow=1,ncol=20)
DD=matrix(0,nrow=1,ncol=20)

for(i in 1:20){
  s.t=s[301,((i-1)*5000+1):((i-1)*5000+5000)]
    E[i]=mean(s.t)
    VA[i]=var(s.t)
    K[i]=mean((s.t-E[i])^4)
    dd=c()
for(j in 1:5000){
  v=s[,((i-1)*5000+j)]
  dd[j]=max(v)-min(v)
}
    DD[i]=mean(dd)
}

DD
ptm<-proc.time()
proc.time()-ptm



#######################################################################################################################
##
## OPtimize cost function 1
##
#######################################################################################################################


library(nloptr)

OPT=function(gamma){

eval_f <- function( x ) {
  g=c()
  for(i in 1:20){
    g[i]=-E[i]+gamma*2*x[i]*VA[i]
  }
  fun=0
  for(i in 1:20){
    fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]
  }
  return( list( "objective" =fun,
                "gradient" = g ) )
}

eval_g_eq <- function( x ) {
  constr <- c( x%*%price -10000 )
  grad <- c( price)
  return( list( "constraints"=constr, "jacobian"=grad ) )
}

# initial values
x0 <- rep(1,20)
# lower and upper bounds of control
lb <- -2500/price[1:20]
ub <- 2500/price[1:20]
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-19 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-19,
              "maxeval" = 10000,
              "local_opts" = local_opts )
res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
return(res)}



OPT(0.5)$solution




#############################################################################################################
##
##  Gamma=0.5
##
##
##
##
##############################################################################################################

POTF=matrix(OPT(0.5)$solution,nrow=20)
library(PerformanceAnalytics)
x <- as.yearmon(2000 + seq(0, 299)/12)
Riskless <- read.csv("~/Desktop/Project/prices/Riskless.csv")
RF=as.numeric(Riskless[,-c(1:301)])/36500
Rf=data.frame(RF,row.names=as.Date(x))
mean(RF)

kl=function(gamma){
POTF=matrix(OPT(gamma)$solution,nrow=20)
N=100
SR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N)
for(i in 1:N){
  ss=c()
  for(j in 1:20){
    ss[j]=(j-1)*5000+round(runif(1,min=1,max=5000),0)
  }
sim=s[,ss]
k=sim%*%POTF
ret=diff(k)/k[1:300] #return
R=data.frame(ret,row.names=as.Date(x))
SR[i]=SharpeRatio(R,Rf,FUN="StdDev") # 
ADD[i]=AverageDrawdown(R) #
MDD[i]=maxDrawdown(R)#
NDD[i]=length(findDrawdowns(R)$from)#
}
return(list(mean(SR),mean(ADD),mean(MDD),mean(NDD)))}



kl(0.0008)


PDF=function(gamma){
POTF=matrix(OPT(gamma)$solution,nrow=20)
TV=numeric(5000)
for(i in 1:5000){
  ss=c()
  for(j in 1:20){
    ss[j]=(j-1)*5000+i
  }
  sim=s[,ss]
  k=sim%*%POTF
  TV[i]=k[301]
}
return(TV)
}


ss=c()
for(j in 1:20){
  ss[j]=(j-1)*5000+5
}




plot(density(PDF(0.5)),xlim=c(9000,30000))
lines(density(PDF(0.1)),col=2)
lines(density(PDF(0.01)),col=3)
lines(density(PDF(0.005)),col=4)
lines(density(PDF(0.003)),col=5)
lines(density(PDF(0.002)),col=6)
lines(density(PDF(0.001)),col=7)
lines(density(PDF(0.0005)),col=8)
lines(density(PDF(0.0001)),col=9)
abline(v=10045.44,lty=2,col=11)


PDF(0.1)
ptm<-proc.time()
proc.time()-ptm

s[,ss]


300*mean(RF)*10000+10000





#######################################################################################################################
##
## OPtimize cost function 1
##
#######################################################################################################################


library(nloptr)

OPT2=function(gamma,lambda){
  
  eval_f <- function( x ) {
    g=c()
    for(i in 1:20){
      g[i]=-E[i]+gamma*2*x[i]*VA[i]+lambda*4*x[i]^3*K[i]
    }
    fun=0
    for(i in 1:20){
      fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]+lambda*x[i]^4*K[i]
    }
    return( list( "objective" =fun,
                  "gradient" = g ) )
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -10000 )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values
  x0 <- rep(1,20)
  # lower and upper bounds of control
  lb <- -2500/price[1:20]
  ub <- 2500/price[1:20]
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-19 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-19,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}


PDF2=function(gamma,lambda){
  POTF=matrix(OPT2(gamma,lambda)$solution,nrow=20)
  TV=numeric(5000)
  for(i in 1:5000){
    ss=c()
    for(j in 1:20){
      ss[j]=(j-1)*5000+i
    }
    sim=s[,ss]
    k=sim%*%POTF
    TV[i]=k[301]
  }
  return(TV)
}

z=0.0000000001
a=rep(0.1,8)
b=c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.1^8)
plot(density(PDF2(0.5,1)),xlim=c(9000,30000))
plot(density(PDF(0.002)),xlim=c(9000,30000))
for(i in 1:8){
  lines(density(PDF2(a[i],b[i])),col=i+1)
}
lines(density(PDF2(0.1,z)),col=2)
lines(density(PDF2(0.01,z)),col=3)
lines(density(PDF2(0.005,z)),col=4)
lines(density(PDF2(0.003,z)),col=5)
lines(density(PDF2(0.002,z)),col=6)
lines(density(PDF2(0.001,z)),col=7)
lines(density(PDF2(0.5,0.1^7)),col=8)
lines(density(PDF2(0.000001,0.1^9*3)),lty=2,col=3)
abline(v=10045.44,lty=2,col=11)

lines(density(PDF2(,1)),col=2)

OPT2(0.02,0)$solution
OPT(0.02)$solution
