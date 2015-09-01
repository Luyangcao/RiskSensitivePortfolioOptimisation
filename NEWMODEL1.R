PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
DATA.T1<-PRICE100[2817:3116,16:19] # Data t-1
DATA.T<-PRICE100[2818:3117,16:19] # Data t

R.T1<-PRICE100[2317:2816,16:19] # Data t-1
R.T<-PRICE100[2318:2817,16:19] # Data t
In.price=as.numeric(DATA.T1[1,])


#######################################################################################################################
##
## fit model X_t=AX_t-1+tau N(0,1)
##
#######################################################################################################################
A=matrix(0,nrow=1,ncol=4)
sd=matrix(0,nrow=1,ncol=4)


for(i in 1:4){
  y=R.T[,i]
  x=R.T1[,i]
  L.Model=lm(y~x-1)
  A[1,i]=as.numeric(L.Model$coefficients)
  sd[1,i]=sd(L.Model$residuals)}

A
sd
plot(L.Model)
plot(density(L.Model$residuals))
#######################################################################################################################
##
## simulations
##
#######################################################################################################################



# 3000* 4 stocks 'true' price
s=matrix(0,nrow=301,ncol=2000*4) # 3min run time
for(i in 1:4){
  for(m in 1:2000){
    s[1,(i-1)*2000+m]=In.price[i] # initial price
    for(n in 1:300){
      s[n+1,(i-1)*2000+m]=A[i]*s[n,(i-1)*2000+m]+rnorm(1,0,sd[i])
    }
  }}



# 
MC.sim=function(IV){
  price.sim=matrix(0,nrow=101,ncol=20000*4)
  for(i in 1:4){
    for(m in 1:20000){
      price.sim[1,(i-1)*20000+m]=IV # initial price
      for(n in 1:100){
        price.sim[n+1,(i-1)*20000+m]=A[i]*price.sim[n,(i-1)*20000+m]+rnorm(1,0,sd[i])
      }
    }}
  E=matrix(0,nrow=1,ncol=4)
  VA=matrix(0,nrow=1,ncol=4)
  K=matrix(0,nrow=1,ncol=4)
  DD=matrix(0,nrow=1,ncol=4)
  for(i in 1:4){
    s.t=price.sim[101,((i-1)*20000+1):((i-1)*20000+20000)]
    E[i]=mean(s.t)
    VA[i]=var(s.t)
    K[i]=mean((s.t-E[i])^4)
    dd=c()
    for(j in 1:20000){
      v=price.sim[,((i-1)*20000+j)]
      dd[j]=max(v)-min(v)
    }
    DD[i]=mean(dd)
  }
  return(list('E'=E,'VA'=VA,'K'=K,'DD'=DD))}

ptm<-proc.time()
MC=MC.sim(1)
proc.time()-ptm
MC
MC2=MC.sim(100)
In.price*A^15

#######################################################################################################################
##
## OPtimize cost function 1
##
#######################################################################################################################


library(nloptr)

OPT1=function(gamma,IH,price){
  E=MC$E*price
  VA=MC$VA
  eval_f <- function( x ) {
    g=c()
    for(i in 1:4){
      g[i]=-E[i]+gamma*2*x[i]*VA[i]
    }
    fun=0
    for(i in 1:4){
      fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]
    }
    return( list( "objective" =fun,
                  "gradient" = g ) )
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  lb <- -(0.75*IH)/abs(price[1:4])
  ub <- (0.75*IH)/abs(price[1:4])
  x0 <- runif(4,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-10 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-10,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}




ptm<-proc.time()
OPT1(0.5,10000,In.price)
proc.time()-ptm

GOPT1<-function(gamma,IH,price){
  set.seed(123)
  K<-OPT1(gamma,IH,price)
  for(i in 1:2){
    set.seed(i+100)
    K1<-OPT1(gamma,IH,price)
    if(K1$objective<K$objective){
      K<-K1
    }
  }
  return(K)
}


fq<-function(x) -5*x^4-3*x^2+6*x





#############################################################################################################
##
##  Gamma=0.5
##
##
##
##
##############################################################################################################

port.sim=function(gamma,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=c()
    for(j in 1:4){
      ss[j]=(j-1)*2000+n
    }
    price=s[,ss]
    
    portfolio=matrix(0,nrow=3,ncol=4)
    portfolio[1,]=OPT1(gamma,10000,In.price)$solution
    X[1:101,n]=price[1:101,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[100*i+1]
      c.price=price[100*i+1,]
      portfolio[i+1,]=OPT1(gamma,IH,c.price)$solution  
      X[(i*100+2):((i+1)*100+1),n]=price[(i*100+2):((i+1)*100+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}



ptm<-proc.time()
X=port.sim(0.5,1)
proc.time()-ptm


pdf.002<-port.sim(0.002,300)




#############################################################################################################
##
##  performance analytics
##
##
##
##############################################################################################################


library(PerformanceAnalytics)
z <- as.yearmon(1960 + seq(0, 299)/12)
Riskless <- read.csv("~/Desktop/Project/prices/Riskless.csv")
RF=as.numeric(Riskless[,-c(1:301)])/36500
Rf=data.frame(RF,row.names=as.Date(z))
mean(RF)

perf=function(gamma,N){
  X=port.sim(gamma,N)
  N=10
  SR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N);MT=numeric(N);AT=numeric(N);Re=numeric(N)
  for(i in 1:N){
    k=p.20.0[,1]
    ret=diff(k)/k[1:300] #return
    R=data.frame(ret,row.names=as.Date(z))
    SR[1]=SharpeRatio(R,Rf,FUN="StdDev") # 
    ADD[i]=AverageDrawdown(R) #
    MDD[i]=maxDrawdown(R)#
    FD=findDrawdowns(R)
    NDD[i]=length(FD$from)
    MT[i]=max(FD$to-FD$from)
    AT[i]=mean(FD$to-FD$from)#
    Re[i]=k[301]/10000-1
  }
  return(list('SR'=mean(SR),'ADD'=mean(ADD),'MDD'=mean(MDD),'NDD'=mean(NDD),
              "MT"=mean(MT),"AT"=mean(AT),"R"=mean(Re),"V"=sd(Re)))}


ptm<-proc.time()
perf(0.5)
proc.time()-ptm

perfYY=function(YYY,N){
  X=YYY
  SR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N);MT=numeric(N);AT=numeric(N);Re=numeric(N)
  for(i in 1:N){
    k=X[,i]
    ret=diff(k)/k[1:300] #return
    R=data.frame(ret,row.names=as.Date(z))
    SR[i]=SharpeRatio(R,Rf,FUN="StdDev") # 
    ADD[i]=AverageDrawdown(R) #
    MDD[i]=maxDrawdown(R)#
    FD=findDrawdowns(R)
    NDD[i]=length(FD$from)
    MT[i]=max(FD$to-FD$from)
    AT[i]=mean(FD$to-FD$from)#
    Re[i]=k[301]/10000-1
  }
  return(list('SR'=mean(SR),'ADD'=mean(ADD),'MDD'=mean(MDD),'NDD'=mean(NDD),
              "MT"=mean(MT),"AT"=mean(AT),"R"=mean(Re),"V"=sd(Re)))}

#############################################################################################################
##
##  PDF
##
##
##
##
##############################################################################################################





PDF=function(gamma){
  X=port.sim(gamma)[301,]
}
fpdf2<-perfYY(pdf2,300)
fpdf1<-perfYY(pdf1,300)
fpdf.5<-perfYY(pdf.5,300)
fpdf.1<-perfYY(pdf.1,300)
fpdf.05<-perfYY(pdf.05,300)
fpdf.01<-perfYY(pdf.01,300)
fpdf.005<-perfYY(pdf.005,300)
fpdf.004<-perfYY(pdf.004,300)
fpdf.003<-perfYY(pdf.003,300)
fpdf.002<-perfYY(pdf.002,300)
fpdf.001<-perfYY(pdf.001,300)
fpdf.0005<-perfYY(pdf.0005,300)
fpdf.0001<-perfYY(pdf.0001,300)


pdf2<-port.sim(2,1000)
pdf1<-port.sim(1,1000)
pdf.5<-port.sim(0.5,1000)
pdf.1<-port.sim(0.1,1000)
pdf.05<-port.sim(0.05,1000)
pdf.01<-port.sim(0.01,1000)
pdf.005<-port.sim(0.005,1000)
pdf.004<-port.sim(0.004,1000)
pdf.003<-port.sim(0.003,1000)
pdf.002<-port.sim(0.002,1000)##################
pdf.001<-port.sim(0.001,1000)
pdf.0005<-port.sim(0.0005,1000)
pdf.0001<-port.sim(0.0001,1000)




plot(density(pdf0.5),xlim=c(10800,25000))
lines(rep(mean(pdf0.5),21),seq(0,2e-04,by=1e-05),lty=6)
lines(density(pdf0.1),col=2)
lines(rep(mean(pdf0.1),21),seq(0,2e-04,by=1e-05),lty=6,col=2)
lines(density(pdf0.01),col=3)
lines(rep(mean(pdf0.01),21),seq(0,2e-04,by=1e-05),lty=6,col=3)
lines(density(pdf0.005),col=4)
lines(rep(mean(pdf0.005),21),seq(0,2e-04,by=1e-05),lty=6,col=4)
lines(density(pdf0.003),col=5)
lines(rep(mean(pdf0.003),21),seq(0,2e-04,by=1e-05),lty=6,col=5)
lines(density(pdf0.002),col=6)
lines(rep(mean(pdf0.002),21),seq(0,2e-04,by=1e-05),lty=6,col=6)
lines(density(pdf0.001),col=7.9999,lwd=1.5)
lines(rep(mean(pdf0.001),21),seq(0,2e-04,by=1e-05),lty=6,col=7.99999)
lines(density(pdf0.0005),col=8.9)
lines(rep(mean(pdf0.0005),21),seq(0,2e-04,by=1e-05),lty=6,col=8.9)
lines(density(pdf0.0001),col=1,lwd=1)
lines(rep(mean(pdf0.0001),21),seq(0,2e-04,by=1e-05),lty=6,col=0,lwd=2)

library(ggplot2)

x1<-c(pdf.5[301,],pdf.1[301,],pdf.01[301,],pdf.005[301,],
      pdf.003[301,],pdf.002[301,],pdf.001[301,],pdf.0005[301,],pdf.0001[301,])
x2=c(rep("0.5",1000),rep("0.1",1000),rep("0.01",1000),rep("0.005",1000),rep("0.003",1000),rep("0.002",1000),rep("0.001",1000),rep("0.0005",1000),rep("0.0001",1000))
k<-data.frame(x1,x2)
colnames(k) <- c("Dollar", "Gamma")
p<-ggplot(data=k,aes(x=Dollar,color=Gamma))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
#p<-p+labs(title='Density of Final Wealth For Different Gamma')
library(plyr)
mu <- ddply(k, "Gamma", summarise, grp.mean=mean(Dollar))
head(mu)
p<-p+geom_vline(data=mu, aes(xintercept=grp.mean, color=Gamma),linetype="dashed")
p






#######################################################################################################################
##
## OPtimize cost function 2
##
#######################################################################################################################


library(nloptr)

OPT2=function(gamma,lambda,IH,price){
  E=MC$E*price
  VA=MC$VA
  K=MC$K
  eval_f <- function( x ) {
    fun=0
    for(i in 1:4){
      p1=0
      for(k in 1:4){
        p1=p1+x[i]^2*VA[i]*x[k]^2*VA[k]
      }
      fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]+lambda*x[i]^4*K[i]+lambda*3*(p1-x[i]^2*VA[i]*x[i]^2*VA[i])
    }
    
    p=0
    for(i in 1:4){
      p=p+x[i]^2*VA[i]
    }
    g=c()
    for(i in 1:4){
      gk=2*x[i]*VA[i]*(p-x[i]^2*VA[i])
      g[i]=-E[i]+gamma*2*x[i]*VA[i]+lambda*4*x[i]^3*K[i]+6*lambda*gk
    }
    
    return( list( "objective" =fun,
                  "gradient" = g ) )
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.75*IH)/abs(price)
  ub <- (0.75*IH)/abs(price)
  set.seed(998)
  x0 <- runif(4,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-19 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-19,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}



OPT2=function(gamma,lambda,IH,price){
  E=MC$E*price
  VA=MC$VA
  K=MC$K
  eval_f <- function( x ) {
    fun=0
    for(i in 1:10){
      p1=0
      for(k in 1:10){
        p1=p1+x[i]^2*VA[i]*x[k]^2*VA[k]
      }
      fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]
    }
    
    p=0
    for(i in 1:10){
      p=p+x[i]^2*VA[i]
    }
    g=c()
    for(i in 1:10){
      gk=2*x[i]*VA[i]*(p-x[i]^2*VA[i])
      g[i]=-E[i]+gamma*2*x[i]*VA[i]
    }
    
    return( list( "objective" =fun,
                  "gradient" = g ) )
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.5*IH)/abs(price)
  ub <- (0.5*IH)/abs(price)
 # set.seed(10000)
  x0 <- runif(10,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-19 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-19,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}


GOPT2<-function(gamma,lambda,IH,price){
  set.seed(123)
K<-OPT2(gamma,lambda,IH,price)
for(i in 1:3){
  set.seed(i+100)
  K1<-OPT2(gamma,lambda,IH,price)
  if(K1$objective<K$objective){
    K<-K1
  }
}
return(K)
}



ptm<-proc.time()
OPT2(0.0001,0,10000,In.price)$solution
proc.time()-ptm

port.sim2=function(gamma,lambda,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=c()
    for(j in 1:4){
      ss[j]=(j-1)*2000+n
    }
    price=s[,ss]
    
    portfolio=matrix(0,nrow=3,ncol=4)
    portfolio[1,]=OPT2(gamma,lambda,10000,In.price)$solution
    X[1:101,n]=price[1:101,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[100*i+1]
      c.price=price[100*i+1,]
      portfolio[i+1,]=OPT2(gamma,lambda,IH,c.price)$solution  
      X[(i*100+2):((i+1)*100+1),n]=price[(i*100+2):((i+1)*100+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}





N=1
ptm<-proc.time()
port.sim2(0.0001,0,1)
proc.time()-ptm





p.5.85<-port.sim2(0.5,5*0.1^8,300)
p.1.85<-port.sim2(0.1,5*0.1^8,300)
p.01.85<-port.sim2(0.01,5*0.1^8,300)
p.001.85<-port.sim2(0.001,5*0.1^8,300)
p.0001.85<-port.sim2(0.0001,5*0.1^8,300)
p.00001.85<-port.sim2(0.00001,5*0.1^8,300)



p.5.6<-p.5.000001
p.1.6<-p.1.000001
p.01.6<-p.01.000001
p.001.6<-p.001.000001
p.0001.6<-p.0001.000001
p.00001.6<-p.00001.000001

plot(density(p.5.6[301,]),xlim=c(10800,25000))

lines(density(p.1.6[301,]),col=2)

lines(density(p.01.6[301,]),col=3)

lines(density(p.001.6[301,]),col=4)

lines(density(p.0001.6[301,]),col=5)

lines(density(p.00001.6[301,]),col=6)


plot(density(p.5.7[301,]),xlim=c(10800,25000))

lines(density(p.1.7[301,]),col=2)

lines(density(p.01.7[301,]),col=3)

lines(density(p.001.7[301,]),col=4)

lines(density(p.0001.7[301,]),col=5)

lines(density(p.00001.7[301,]),col=6)



plot(density(p.5.85[301,]),xlim=c(10800,25000))

lines(density(p.1.85[301,]),col=2)

lines(density(p.01.85[301,]),col=3)

lines(density(p.001.85[301,]),col=4)

lines(density(p.0001.85[301,]),col=5)

lines(density(p.00001.85[301,]),col=6)


plot(density(p.5.8[301,]),xlim=c(10800,25000))

lines(density(p.1.8[301,]),col=2)

lines(density(p.01.8[301,]),col=3)

lines(density(p.001.8[301,]),col=4)

lines(density(p.0001.8[301,]),col=5)

lines(density(p.00001.8[301,]),col=6)


plot(density(p.5.11[301,]),xlim=c(10800,25000))

lines(density(p.1.11[301,]),col=2)

lines(density(p.01.11[301,]),col=3)

lines(density(p.001.11[301,]),col=4)

lines(density(p.0001.11[301,]),col=5)

lines(density(p.00001.11[301,]),col=6)

OPT2(0.5,0.1^9,10000,In.price)


p.002.1<-port.sim2(0.002,0.1,300)
p.002.2<-port.sim2(0.002,0.1^2,300)
p.002.3<-port.sim2(0.002,0.1^3,300)
p.002.4<-port.sim2(0.002,0.1^4,300)
p.002.5<-port.sim2(0.002,0.1^5,300)
p.002.6<-port.sim2(0.002,0.1^6,300)
p.002.7<-port.sim2(0.002,0.1^7,300)
p.002.8<-port.sim2(0.002,0.1^8,300)#
p.002.9<-port.sim2(0.002,0.1^9,300)
p.002.10<-port.sim2(0.002,0.1^10,300)




plot(density(p.5.1[301,]),xlim=c(10800,25000))

lines(density(p.5.01[301,]),col=2)

lines(density(p.5.0001[301,]),col=3)

lines(density(p.5.000001[301,]),col=4)

lines(density(p.5.7[301,]),col=5)

lines(density(p.5.85[301,]),col=6)
lines(density(p.5.9[301,]),col=1)
lines(density(p.5.11[301,]),col=2)
lines(density(p.5.15[301,]),col=1)
lines(density(p.5.19[301,]),col=1)





plot(density(p.1.1[301,]),xlim=c(10800,25000))

lines(density(p.1.01[301,]),col=2)

lines(density(p.1.0001[301,]),col=3)

lines(density(p.1.000001[301,]),col=4)

lines(density(p.1.7[301,]),col=5)

lines(density(p.1.85[301,]),col=6)
lines(density(p.1.8[301,]),col=6)
lines(density(p.1.9[301,]),col=1)
lines(density(p.1.11[301,]),col=2)
lines(density(p.1.15[301,]),col=1)
lines(density(p.1.19[301,]),col=1)



plot(density(p.01.1[301,]),xlim=c(10800,25000))
lines(density(p.01.01[301,]),col=2)

lines(density(p.01.0001[301,]),col=3)

lines(density(p.01.000001[301,]),col=4)

lines(density(p.01.7[301,]),col=5)

lines(density(p.01.85[301,]),col=6)
lines(density(p.01.8[301,]),col=6)
lines(density(p.01.9[301,]),col=1)
lines(density(p.01.11[301,]),col=2)
lines(density(p.01.15[301,]),col=1)
lines(density(p.01.19[301,]),col=1)


plot(density(p.001.1[301,]),xlim=c(10000,20000))
lines(density(p.001.6[301,]),col=2)
lines(density(p.001.7[301,]),col=3)
lines(density(p.001.7[301,]),col=3)
lines(density(p.001.8[301,]),xlim=c(10000,20000))
lines(density(p.001.95[301,]),col=2)
lines(density(p.001.9[301,]),col=3)
lines(density(p.001.10[301,]),col=3)
lines(density(p.001.11[301,]),col=4)
lines(density(p.001.12[301,]),col=5)
lines(density(p.001.13[301,]),col=6)
lines(density(p.001.14[301,]),col=7)
lines(density(p.001.15[301,]),col=8)
lines(density(p.001.18[301,]),col=9)
lines(density(p.001.20[301,]),col=9)
lines(density(p.001.0[301,]),col=9)
