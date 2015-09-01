library(PerformanceAnalytics)
library(nloptr)
library(Rcpp)
library(inline)
library(ggplot2)




ptm<-proc.time()
port.simM25(6000,2)
proc.time()-ptm


########################################################################################################################

DATASET=function(){
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
  DATA.T1<-PRICE100[2817:3116,2:51] # Data t-1
  DATA.T<-PRICE100[2818:3117,2:51] # Data t
  
  R.T1<-PRICE100[2317:2816,2:51] # Data t-1
  R.T<-PRICE100[2318:2817,2:51] # Data t
  In.price=as.numeric(DATA.T1[1,])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
}

########################################################################################################################

Fit_Model=function(){
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  B=matrix(0,nrow=1,ncol=50)
  sig=matrix(0,nrow=1,ncol=50)
  
  
  for(i in 1:50){
    y=R.T[,i]
    x=R.T1[,i]
    L.Model=lm(log(y/x)~1)
    B[1,i]=as.numeric(L.Model$coefficients)
    sig[1,i]=sd(L.Model$residuals)}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
}

########################################################################################################################

Sim_TruePath=function(){
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cxxsm2<-'
Environment stats("package:stats");
  Function rnorm = stats["rnorm"];
NumericMatrix sm2(301,50000);

  NumericVector B(RB);
  NumericVector sig(Rsig);
  NumericVector  price(Rprice);
  
  for (int i=1; i<=50; i++) {
  for(int m=1; m<=1000; m++){
  sm2(1-1,(i-1)*1000+m-1)=price[i-1];
  for(int n=1; n<=300; n++){
  sm2(n+1-1,(i-1)*1000+m-1)=sm2(n-1,(i-1)*1000+m-1)*exp(B(i-1)+as<double>(rnorm(1,0,sig[i-1])));
  }
  }
  }
  return wrap(sm2);'
  
  Cppsm2<-cxxfunction(signature(RB='vector',Rsig='vector',Rprice='vector'),
                    body=cxxsm2, plugin="Rcpp")
  sm2=Cppsm2(B,sig,In.price)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
}

########################################################################################################################

MC_sim=function(){
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Mcsimbodym2<-'
  Environment stats("package:stats");
  Environment base("package:base");
  Function mean = base["mean"];
  Function var = stats["var"];
  Function rnorm = stats["rnorm"];
  Function max = base["max"];
  Function min = base["min"];

  NumericVector B(RB);
  NumericVector sig(Rsig);
  
  double IV = as<double>(RIV);
  NumericMatrix pricesim(61,10000*50);
  
  
  
  for(int i=1;i<=50;i++){
  for(int m=1; m<=10000;m++){
  pricesim(1-1,(i-1)*10000+m-1)=IV;
  for(int n=1;n<=60;n++){
  pricesim(n+1-1,(i-1)*10000+m-1)=pricesim(n-1,(i-1)*10000+m-1)*exp(B[i-1]+as<double>(rnorm(1,0,sig[i-1])));
  }
  }}
  
  NumericVector E(50);
  NumericVector VA(50);
  NumericVector K(50);
  NumericVector DD(50);
  NumericVector v(61);

  for(int i=1; i<=50;i++){
  NumericVector st(10000);
  NumericVector dd(10000);
  for(int m=1; m<=10000;m++){
  st[m-1]=pricesim(61-1,(i-1)*10000+m-1);
  v=pricesim(_,(i-1)*10000+m-1);
  dd[m-1]=as<double>(max(v))-as<double>(min(v));
  }
  DD[i-1]=as<double>(mean(dd));
  E[i-1]=as<double>(mean(st));
  VA[i-1]=as<double>(var(st));
  NumericVector st2(10000);
  for(int n=1; n<=10000;n++){
  st2[n-1]=(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1]);}
  K[i-1]=as<double>(mean(st2));
  }
  List ret; ret["E"] = E; ret["VA"] = VA;ret["K"] = K;ret["DD"] = DD;
  return wrap(ret);
  '
  wwm2<-cxxfunction(signature(RIV='numeric',RB='vector',Rsig='vector'),
                  body=Mcsimbodym2, plugin="Rcpp")
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MCm2<-wwm2(1,B,sig)
  MCm2.100<-wwm2(100,B,sig)
  
  
}

########################################################################################################################
%%%
EVALFF<-function(){
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  evalbody4<-'
  double gamma = as<double>(Rgamma);
  double lambda = as<double>(Rlambda);
  double theta = as<double>(Rtheta);
  NumericVector x(Rx);
  NumericVector E(RE);
  NumericVector VA(RVA);
  NumericVector K(RK);
  NumericVector DD(RDD);
  
  NumericVector fun(1);
  NumericVector g(50);
  NumericVector p2(1);
  NumericVector gk(1);
  
  for(int i=1; i<=50; i++){
  NumericVector p1(1);
  for(int k=1; k<=50; k++){
  p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1);
  }
  fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)+lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1))+theta*x(i-1)*DD(i-1);
  }
  
  for(int i=1;i<=50;i++){
  p2=p2+x(i-1)*x(i-1)*VA(i-1);
  }
  
  
  for(int i=1; i<=50; i++){
  gk=2*x(i-1)*VA(i-1)*(p2-x(i-1)*x(i-1)*VA(i-1));
  g(i-1)=-E(i-1)+gamma*2*x(i-1)*VA(i-1)+lambda*4*x(i-1)*x(i-1)*x(i-1)*K(i-1)+6*lambda*gk(0)+theta*DD(i-1);
  }
  
  
  List ret; ret["objective"] = fun; ret["gradient"] = g;
  return wrap(ret);
  '
  eval_f4<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rtheta='numeric',Rx='vector',RE='vector',RVA='vector',RK='vector',RDD='vector'),
                       body=evalbody4, plugin="Rcpp")
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
}


OPTM2=function(gamma,lambda,theta,IH,price){
  E=MCm2$E*price
  VA=MCm2$VA*(price^2)
  K=MCm2$K*(price^4)
  DD=MCm2$DD*price
  
  eval_f <- function( x ) {eval_f4(gamma,lambda,theta,x,E,VA,K,DD)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.1*IH)/abs(price)
  ub <- (0.1*IH)/abs(price)
  set.seed(998)
  x0 <- runif(50,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-8 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-8,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}

########################################################################################################################
$$$$$
  Cssf<-function(){
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ssbody<-'
    double n = as<double>(Rn);
    NumericVector ss(50);
    for(int i=1; i<=50; i++){
    ss(i-1)=(i-1)*1000+n;
    }
    return wrap(ss);
    '
    Css<-cxxfunction(signature(Rn='numeric'),
                     body=ssbody, plugin="Rcpp")
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  }


port.simM2=function(gamma,lambda,theta,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=sm2[,ss]
    
    portfolio=matrix(0,nrow=5,ncol=50)
    portfolio[1,]=OPTM2(gamma,lambda,theta,10000,In.price)$solution
    X[1:61,n]=price[1:61,]%*%portfolio[1,]
    for(i in 1:4){
      IH=X[60*i+1,n]
      c.price=price[60*i+1,]
      portfolio[i+1,]=OPTM2(gamma,lambda,theta,IH,c.price)$solution  
      X[(i*60+2):((i+1)*60+1),n]=price[(i*60+2):((i+1)*60+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

########################################################################################################################

PerfDATA<-function(){
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  z <- as.yearmon(1960 + seq(0, 299)/12)
  RISKFREE <- read.csv("~/Desktop/Project/prices/RISKFREE.csv")
  RF=as.numeric(RISKFREE[1:300,2])/36500
  Rf=data.frame(RF,row.names=as.Date(z))
  mean(RF)
  
  
  Nasdaq <- read.csv("~/Desktop/Project/prices/Nasdaq.csv")
  bench<-rev(Nasdaq[10:310,c(7)])
  Rben=diff(bench)/bench[1:300]
  Rb=data.frame(Rben,row.names=as.Date(z))
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
}

perfYY=function(YYY,N){
  X=YYY
  SR=numeric(N);IR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N);MT=numeric(N);AT=numeric(N);Re=numeric(N)
  for(i in 1:N){
    k=X[,i]
    ret=diff(k)/k[1:300] #return
    R=data.frame(ret,row.names=as.Date(z))
    SR[i]=SharpeRatio(R,Rb,FUN="StdDev") # 
    ADD[i]=AverageDrawdown(R) #
    MDD[i]=maxDrawdown(R)#
    FD=findDrawdowns(R)
    NDD[i]=length(FD$from)
    MT[i]=max(FD$to-FD$from)
    AT[i]=mean(FD$to-FD$from)#
    Re[i]=k[301]/10000-1
    IR[i]=InformationRatio(R,Rb)
  }
  return(list('SR'=mean(SR),'IR'=mean(IR),'MDD'=mean(MDD),'ADD'=mean(ADD),'NDD'=mean(NDD),
              "MT"=mean(MT),"AT"=mean(AT),"R"=mean(Re),"V"=sd(Re)))}


########################################################################################################################

########################################################################################################################
