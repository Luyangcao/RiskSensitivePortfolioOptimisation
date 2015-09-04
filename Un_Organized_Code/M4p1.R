res=list() #4,9,15,18,28,29,31,46
BBB=BBlack[,-c(4,9,15,18,23,28,29,31,41,46)]
for(i in 21:40){
  y=BBB[1:100,i]
  rett=log(y+1)
  res[[i-20]]=ugarchfit(data =rett, spec = ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE)))
  }


library(rugarch)

MKM=function(){
SIM=matrix(0,nrow=480,ncol=20000)
for(i in 1:20){
sim=ugarchsim(res[[i]], n.sim = 480, m.sim = 1000, n.start=100)
JK=matrix(0,nrow=480,ncol=1000)
for(j in 1:1000){
JK[,j]=cumsum(fitted(sim)[,j])}
SIM[,((i-1)*1000+1):(i*1000)]=exp(JK)
}
return(SIM)
}
HAHA=MKM()
HAHA=rbind(rep(1,20000),HAHA)




SSIM=function(){
  SIM=matrix(0,nrow=12,ncol=20000)
  for(i in 1:20){
    sim=ugarchsim(res[[i]], n.sim = 12, m.sim = 1000, n.start=100)
    JK=matrix(0,nrow=12,ncol=1000)
    for(j in 1:1000){
      JK[,j]=cumsum(fitted(sim)[,j])}
    SIM[,((i-1)*1000+1):(i*1000)]=exp(JK)
  }
  return(SIM)
}
SIM=SSIM()
DF=SSIM()
DF1=SSIM()
DF2=SSIM()
DF3=SSIM()
DF4=SSIM()
DF5=SSIM()
?ugarchsim

property=function(seed,N){
MCE=c();MCVA=c();MCK=c()
for(i in 1:20){
  sim=ugarchsim(res[[i]], n.sim = 12, m.sim = N, n.start=100,rseed = seed)
  JK=colSums(fitted(sim))
  M=mean(exp(JK))
  while(M>1.2){
    sim=ugarchsim(res[[i]], n.sim = 12, m.sim = N, n.start=100,rseed = seed)
    JK=colSums(fitted(sim))
    M=mean(exp(JK))
  }
  MCE[i]=mean(exp(JK))
  MCVA[i]=var(exp(JK))
  MCK[i]=mean((exp(JK)-MCE[i])^4)
}
return(list('E'=MCE,'VA'=MCVA,'K'=MCK))
}
MC=property(NA,10000)
MC1=property(NA,10000)
MC2=property(NA,10000)
MC3=property(NA,10000)
MC2=MC
MC=MC1
pwww=function(){
MCE=c();MCVA=c();MCK=c()
for(i in 1:40){
KL=(SIM[12,((i-1)*10000+1):(i*10000)])
MCE[i]=mean((KL))
MCVA[i]=var((KL))
MCK[i]=mean(((KL)-MCE[i])^4)
}
return(list('E'=MCE,'VA'=MCVA,'K'=MCK))
}



MME=c();MMVA=c()
for(i in 1:40){
ff=ugarchforecast(res[[i]],  n.ahead = 36)
MME[i]=exp(cumsum(fitted(ff)))[12]
MMVA[i]=cumsum(sigma(ff)^2)[12]}

MMVA
MC$VA

dim(fitted(sim))
MC=MC2
MC=colSums(fitted(sim))
length(MC)
ptm<-proc.time()
MC=property(1:100000,100000)
proc.time()-ptm


evalbody4<-'
double gamma = as<double>(Rgamma);
double lambda = as<double>(Rlambda);

NumericVector x(Rx);
NumericVector E(RE);
NumericVector VA(RVA);
NumericVector K(RK);


NumericVector fun(1);
NumericVector g(20);
NumericVector p2(1);
NumericVector gk(1);

for(int i=1; i<=20; i++){
NumericVector p1(1);
for(int k=1; k<=20; k++){
p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1);
}
fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)+lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1));
}

for(int i=1;i<=20;i++){
p2=p2+x(i-1)*x(i-1)*VA(i-1);
}


for(int i=1; i<=20; i++){
gk=2*x(i-1)*VA(i-1)*(p2-x(i-1)*x(i-1)*VA(i-1));
g(i-1)=-E(i-1)+gamma*2*x(i-1)*VA(i-1)+lambda*4*x(i-1)*x(i-1)*x(i-1)*K(i-1)+6*lambda*gk(0);
}


List ret; ret["objective"] = fun; ret["gradient"] = g;
return wrap(ret);
'
eval_f4<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rx='vector',RE='vector',RVA='vector',RK='vector'),
                     body=evalbody4, plugin="Rcpp")





OPTM3=function(gamma,lambda,IH,price){
  
  E=MC$E*price
  VA=MC$VA*price^2
  K=MC$K*price^4
  
  eval_f <- function( x ) {eval_f4(gamma,lambda,x,E,VA,K)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.3*IH)/abs(price)
  ub <- (0.3*IH)/abs(price)
  set.seed(998)
  x0 <- runif(20,lb,0.5/3*lb)
  x1 <- runif(20,0.5/3*lb,0.25/3*ub)
  x2 <- runif(20,0.25/3*ub,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-8 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-8,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res0 <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_ineq=eval_g_eq,opts=opts)
  res1 <- nloptr(x1,eval_f=eval_f,lb=lb,ub=ub,eval_g_ineq=eval_g_eq,opts=opts)
  res2 <- nloptr(x2,eval_f=eval_f,lb=lb,ub=ub,eval_g_ineq=eval_g_eq,opts=opts)
  if(res0$objective<=res1$objective){
    res=res0
  }else{res=res1}
  if(res$objective<=res2$objective){
    res=res
  }else{res=res2}
  return(res)}



ssbody<-'
    double n = as<double>(Rn);
NumericVector ss(20);
for(int i=1; i<=20; i++){
ss(i-1)=(i-1)*1000+n;
}
return wrap(ss);
'
Css<-cxxfunction(signature(Rn='numeric'),
                 body=ssbody, plugin="Rcpp")




port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=37,ncol=N)
  X[1,]=rep(1,N)
  portfolio=matrix(0,nrow=3,ncol=40)
  portfolio[1,]=OPTM3(gamma,lambda,1,rep(1,40))$solution
  weight=portfolio[1,]*rep(1,40)/1
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=SIM[,ss]
    
    X[2:13,n]=price[1:12,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[12*i+1,n]
      c.price=price[12*i,]
      portfolio[i+1,]=weight*IH/c.price
      X[(i*12+2):((i+1)*12+1),n]=price[(i*12+1):((i+1)*12),]%*%portfolio[i+1,]
    }
  }
  return(X)
}
