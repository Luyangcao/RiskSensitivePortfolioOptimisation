

MC.sim3=function(price,N){
  price=price
  price.sim=matrix(0,nrow=61,ncol=N*50)
  for(i in 1:50){
    for(m in 1:N){
      price.sim[1,(i-1)*N+m]=price[i] # initial price
      for(n in 1:60){
        price.sim[n+1,(i-1)*N+m]=price.sim[n,(i-1)*N+m]*exp(B[i]+rnorm(1,0,sig[i]))
      }
    }}

  DD=matrix(0,nrow=1,ncol=50)
  for(i in 1:50){
    dd=c()
    for(j in 1:N){
      v=price.sim[,((i-1)*N+j)]
      dd[j]=max(v)-min(v)
    }
    DD[i]=mean(dd)
  }
  return(list('DD'=DD))}


MC.sim3(rep(1000,50),100)$DD






MCDD1=as.numeric(MC.sim3(In.price,30000)$DD)

MCDD2=matrix(0,nrow=300,ncol=4)
MCDD3=matrix(0,nrow=300,ncol=4)
for(n in 1:300){
  ss=c()
  for(j in 1:4){
    ss[j]=(j-1)*2000+n
  }
  price2=s[101,ss]
  price3=s[201,ss]
  MCDD2[n,]=as.numeric(MC.sim3(price2,2000)$DD)
  MCDD3[n,]=as.numeric(MC.sim3(price3,2000)$DD)
  }
ptm<-proc.time()
MC.sim3(In.price,2000)$DD
proc.time()-ptm


OPT3=function(gamma,theta,IH,price,DrawD){
  E=MC$E*price
  VA=MC$VA
  DD=as.numeric(DrawD)
  eval_f <- function( x ) {
    g=c()
    for(i in 1:4){
      g[i]=-E[i]+gamma*2*x[i]*VA[i]+theta*D[i]
    }
    fun=0
    for(i in 1:4){
      fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]+theta*x[i]*D[i]
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
  lb <- -(0.75*IH)/abs(price[1:4])
  ub <- (0.75*IH)/abs(price[1:4])
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

OPT3(0.5,0.1,10000,In.price,MCDD1)

D=as.numeric(MCDD1[1,])
D[8]

GOPT3<-function(gamma,theta,IH,price,DD){
  K<-OPT3(gamma,theta,IH,price,DD)
   for(i in 1:2){
     K1<-OPT2(gamma,theta,IH,price)
     if(K1$objective<K$objective){
      K<-K1
     }
 }
  return(K)
}



port.sim3=function(gamma,theta,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=c()
    for(j in 1:4){
      ss[j]=(j-1)*2000+n
    }
    price=s[,ss]
    
    portfolio=matrix(0,nrow=3,ncol=4)
      portfolio[1,]=OPT3(gamma,theta,10000,In.price,MCDD1)$solution
      X[1:101,n]=price[1:101,]%*%portfolio[1,]
    
      IH=X[100*1+1]
      c.price=price[100*1+1,]
      portfolio[1+1,]=OPT3(gamma,theta,IH,c.price,MCDD2[n])$solution  
      X[(1*100+2):((1+1)*100+1),n]=price[(1*100+2):((1+1)*100+1),]%*%portfolio[1+1,]
      
      IH=X[100*2+1]
      c.price=price[100*2+1,]
      portfolio[2+1,]=OPT3(gamma,theta,IH,c.price,MCDD2[n])$solution  
      X[(2*100+2):((2+1)*100+1),n]=price[(2*100+2):((2+1)*100+1),]%*%portfolio[2+1,]
  }
  return(X)
}
#######

G0.0006=function(){
  p3.0006.0<-port.sim3(0.0006,0,300)
  p3.0006.20<-port.sim3(0.0006,2,300)
  p3.0006.1<-port.sim3(0.0006,0.1,300)
  p3.0006.2<-port.sim3(0.0006,0.1^2,300)
  p3.0006.3<-port.sim3(0.0006,0.1^3,300)
  p3.0006.4<-port.sim3(0.0006,0.1^4,300)
  p3.0006.5<-port.sim3(0.0006,0.1^5,300)###############
  p3.0006.6<-port.sim3(0.0006,0.1^6,300)
  p3.0006.7<-port.sim2(0.0006,0.1^7,300)
  p3.0006.8<-port.sim2(0.0006,0.1^8,300)#
  p.0006.95<-port.sim2(0.0006,5*0.1^9,300)
  p.0006.9<-port.sim2(0.0006,0.1^9,300)
  p.0006.105<-port.sim2(0.0006,5*0.1^10,300)
  p.0006.10<-port.sim2(0.0006,0.1^10,300)
  p.0006.115<-port.sim2(0.0006,5*0.1^11,300)
  p.0006.11<-port.sim2(0.0006,0.1^11,300)
  p.0006.125<-port.sim2(0.0006,5*0.1^12,300)
  p.0006.12<-port.sim2(0.0006,0.1^12,300)
  p.0006.13<-port.sim2(0.0006,0.1^13,300)
  p.0006.14<-port.sim2(0.0006,0.1^14,300)
  p.0006.15<-port.sim2(0.0006,0.1^15,300)
  p.0006.18<-port.sim2(0.0006,0.1^18,300)
  p.0006.22<-port.sim2(0.0006,0.1^22,300)
  
  
  
  f3.0006.0<-perfYY(p3.0006.0,300)
  f3.0006.20<-perfYY(p3.0006.20,300)
  f3.0006.1<-perfYY(p3.0006.1,300)
  f3.0006.2<-perfYY(p3.0006.2,300)
  f3.0006.3<-perfYY(p3.0006.3,300)
  f3.0006.4<-perfYY(p3.0006.4,300)
  f3.0006.5<-perfYY(p3.0006.5,300)
  f3.0006.6<-perfYY(p3.0006.6,300)
  f.0006.7<-perfYY(p.0006.7,300)
  f.0006.8<-perfYY(p.0006.8,300)
  f.0006.95<-perfYY(p.0006.95,300)
  f.0006.9<-perfYY(p.0006.9,300)
  f.0006.105<-perfYY(p.0006.105,300)
  f.0006.10<-perfYY(p.0006.10,300)
  f.0006.115<-perfYY(p.0006.115,300)
  f.0006.11<-perfYY(p.0006.11,300)
  f.0006.125<-perfYY(p.0006.125,300)
  f.0006.12<-perfYY(p.0006.12,300)
  f.0006.13<-perfYY(p.0006.13,300)
  f.0006.14<-perfYY(p.0006.14,300)
  f.0006.15<-perfYY(p.0006.15,300)
  f.0006.18<-perfYY(p.0006.18,300)
  f.0006.22<-perfYY(p.0006.22,300)
  f3.0006.0
}

rbind(f3.0006.20,f3.0006.1,f3.0006.2,f3.0006.3,f3.0006.4,f3.0006.5,f3.0006.6,f3.0006.0)
plot(density(p3.0006.0[301,]),xlim=c(5000,30000))
lines(density(p3.0006.20[301,]),col=2)
lines(density(p3.0006.1[301,]),col=3)
lines(density(p3.0006.2[301,]),col=3)
lines(density(p3.0006.3[301,]),col=3)






psim1=matrix(0,nrow=10,ncol=10000)
psim1[1,]=5

for(i in 1:10000){
  for(k in 1:9){
    psim1[k+1,i]=1.01*psim1[k,i]+rnorm(1,0,0.5)
  }
}

dddd=c()
for(i in 1:10000){
  k=psim1[,i]
  v=cumsum(k)/c(1:10)
  dddd[i]=max(v)-min(v)
}
mean(dddd)



psim2=matrix(0,nrow=10,ncol=10000)
psim2[1,]=500

for(i in 1:10000){
  for(k in 1:9){
    psim2[k+1,i]=1.01*psim2[k,i]+rnorm(1,0,2)
  }
}
