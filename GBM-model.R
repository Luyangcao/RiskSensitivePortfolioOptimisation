library(PerformanceAnalytics)
library(nloptr)
library(Rcpp)
library(inline)
library(ggplot2)
library(ghyp)
library(rugarch)



########################################################################################################################


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
  DATA.T1<-PRICE100[3027:3116,21:40] # Data t-1
  DATA.T<-PRICE100[3028:3117,21:40] # Data t
  
  R.T1<-PRICE100[2317:3026,21:40] # Data t-1
  R.T<-PRICE100[2318:3027,21:40] # Data t
  In.price=as.numeric(DATA.T1[1,])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

########################################################################################################################


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  B=matrix(0,nrow=1,ncol=20)     # mu-sigma^2/2
  sig=matrix(0,nrow=1,ncol=20)   # sigma
  
  
  for(i in 1:20){
    y=R.T[,i]
    x=R.T1[,i]
    L.Model=lm(log(y/x)~1)    # linear regression
    B[1,i]=as.numeric(L.Model$coefficients)
    sig[1,i]=sd(L.Model$residuals)}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  


########################################################################################################################


  #following code is used to simulate price from GBM model.C++ code

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  price_simulate_body<-'
Environment stats("package:stats"); 
  Function rnorm = stats["rnorm"];
NumericMatrix sm2(91,200000);

  NumericVector B(RB);
  NumericVector sig(Rsig);
  NumericVector  price(Rprice);
  
  for (int i=1; i<=20; i++) {
  for(int m=1; m<=10000; m++){
  sm2(1-1,(i-1)*10000+m-1)=price[i-1];
  for(int n=1; n<=90; n++){
  sm2(n+1-1,(i-1)*10000+m-1)=sm2(n-1,(i-1)*10000+m-1)*exp(B(i-1)+as<double>(rnorm(1,0,sig[i-1])));
  }
  }
  }
  return wrap(sm2);'
  
  PriceSim<-cxxfunction(signature(RB='vector',Rsig='vector',Rprice='vector'),
                    body=price_simulate_body, plugin="Rcpp")

  SimPrice=PriceSim(B,sig,In.price)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
}

########################################################################################################################


  
  #Following code is used to find Mean Vairance and kurtosis for each stock at end of the period. C++ code
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MCsimbody<-'
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
  NumericMatrix pricesim(31,10000*20);
  
  
  
  for(int i=1;i<=20;i++){
  for(int m=1; m<=10000;m++){
  pricesim(1-1,(i-1)*10000+m-1)=IV;
  for(int n=1;n<=30;n++){
  pricesim(n+1-1,(i-1)*10000+m-1)=pricesim(n-1,(i-1)*10000+m-1)*exp(B[i-1]+as<double>(rnorm(1,0,sig[i-1])));
  }
  }}
  
  NumericVector E(20);
  NumericVector VA(20);
  NumericVector K(20);
  NumericVector v(31);

  for(int i=1; i<=20;i++){
  NumericVector st(10000);
  for(int m=1; m<=10000;m++){
  st[m-1]=pricesim(31-1,(i-1)*10000+m-1);
  }
  E[i-1]=as<double>(mean(st));
  VA[i-1]=as<double>(var(st));
  NumericVector st2(10000);
  for(int n=1; n<=10000;n++){
  st2[n-1]=(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1]);}
  K[i-1]=as<double>(mean(st2));
  }
  List ret; ret["E"] = E; ret["VA"] = VA;ret["K"] = K;
  return wrap(ret);
  '
  MCproperty<-cxxfunction(signature(RIV='numeric',RB='vector',Rsig='vector'),
                  body=MCsimbody, plugin="Rcpp")
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MCm2<-MCproperty(1,B,sig)
  MCm2.100<-wwm2(100,B,sig)
  MCm2
  


########################################################################################################################

  #following is reward function mean-variance-kurtosis and its gradient. This function is used in optimisation. 
  
  #C++ code
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  # optimization function for geometric brwonian motion model mean variance kurtosis reward function
  # gamma is parameter of variance, lambda is parameter of kurtosis
  #set lambda=0, the reward function is mean-vairance
  # IH= initial holding,

OPTM2=function(gamma,lambda,IH,price){    
  E=MCm2$E*price
  VA=MCm2$VA*(price^2)
  K=MCm2$K*(price^4)
  
  eval_f <- function( x ) {eval_f4(gamma,lambda,x,E,VA,K)}
  
  eval_g_eq <- function( x ) {    # equality constriant
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.25*IH)/abs(price)
  ub <- (0.25*IH)/abs(price)
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

########################################################################################################################

  
#  a simple function calculate the number of columns, used to confirm column of each stock!
  Cssf<-function(){
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ssbody<-'
    double n = as<double>(Rn);
    NumericVector ss(20);
    for(int i=1; i<=20; i++){
    ss(i-1)=(i-1)*10000+n;
    }
    return wrap(ss);
    '
    Css<-cxxfunction(signature(Rn='numeric'),
                     body=ssbody, plugin="Rcpp")
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  }

#### calculate portfolio value of GBM model of mean-variance kurtoisis reward function
MVK_portfolio=function(gamma,lambda,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=20)
  portfolio[1,]=OPTM2(gamma,lambda,1,In.price)$solution
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=SimPrice[,ss]
  
    X[1:31,n]=price[1:31,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[30*i+1,n] 
      c.price=price[30*i+1,]
      portfolio[i+1,]=OPTM2(gamma,lambda,IH,c.price)$solution  
      X[(i*30+2):((i+1)*30+1),n]=price[(i*30+2):((i+1)*30+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

########################################################################################################################


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  z <- as.yearmon(1960 + seq(0, 299)/12)
 #calculate risk free rate
  RISKFREE <- read.csv("~/Desktop/Project/prices/RISKFREE.csv")
  RF=as.numeric(RISKFREE[1:300,2])/36500
  Rf=data.frame(RF,row.names=as.Date(z))
  mean(RF)
  
  
  # calculate bench of nasdaq market
  Nasdaq <- read.csv("~/Desktop/Project/prices/Nasdaq.csv")
  bench<-rev(Nasdaq[10:310,c(7)])
  Rben=diff(bench)/bench[1:300]
  Rb=data.frame(Rben,row.names=as.Date(z))

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
# calcualte the perfomance indices, SR=sharpe ratio, IR=informantion ratio, ADD=average drawdown, MDD=maximum drawdown
# NDD=number of drawdowan, MT=maximum time of drawdown, LS=loss probability, Re=return, V=volatility

performance=function(YYY,N){
  X=YYY   #portfolio values
  SR=numeric(N);IR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N);MT=numeric(N);LS=numeric(N);Re=numeric(N)
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
    LS[i]=as.numeric(k[301]<1.081733)
    Re[i]=k[301]-1
    IR[i]=InformationRatio(R,Rb)
  }
  return(list('SR'=mean(SR),'IR'=mean(IR),'MDD'=mean(MDD),'ADD'=mean(ADD),'NDD'=mean(NDD),
              "MT"=mean(MT),"LS"=mean(LS),"R"=mean(Re),"V"=sd(Re)))}


########################################################################################################################

########################################################################################################################

  
  
  ###risk sensitive reward function and its gradient, C++ code
  
  ## mcs or mcs1 is simulated terminal value of each stock, used to monte carlo method in optimization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  RS_GBM_body<-'
  #include <math.h>
  
  double xi = as<double>(Rxi);
  
  
  NumericVector x(Rx);
  NumericVector price(Rprice);
  NumericMatrix mcs(Rmcs1);
  
  
  NumericVector fun(1);
  NumericVector g(50);
  NumericVector Vmcs(50);
  NumericVector XP(100);
  
  for(int j=1; j<=100; j++){
  Vmcs=mcs(_,j-1);
  for(int i=1; i<=50; i++){
  XP(j-1)=XP(j-1)+x(i-1)*Vmcs(i-1)*price(i-1);
  }
  fun=fun+exp(-xi*XP(j-1)+xi*10000);
  }
  
  
  
  for(int k=1; k<=50; k++){
  NumericVector gfz(1);
  for(int j=1; j<=100; j++){
  Vmcs=mcs(_,j-1);
  gfz=gfz+exp(-xi*XP(j-1)+xi*10000)*(-xi*Vmcs[k-1]*price(k-1));}
  g(k-1)=as<double>(gfz);
  }
  
  
  List ret; ret["objective"] = fun; ret["gradient"] = g;
  return wrap(ret);
  '
  RS_GBM<-cxxfunction(signature(Rxi='numeric',Rx='vector',Rprice='vector',Rmcs1='Matrix'),
                       body=RS_GBM_body, plugin="Rcpp")
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  
  
  
  
  ## mcs or mcs1 is simulated terminal value of each stock, used to monte carlo method in optimization
  
  mcbody<-'
  Environment stats("package:stats");
  Function rnorm = stats["rnorm"];
  NumericMatrix sm2(61,50000);
  
  NumericVector B(RB);
  NumericVector sig(Rsig);
  NumericVector  price(Rprice);
  
  for (int i=1; i<=50; i++) {
  for(int m=1; m<=1000; m++){
  sm2(1-1,(i-1)*1000+m-1)=price[i-1];
  for(int n=1; n<=60; n++){
  sm2(n+1-1,(i-1)*1000+m-1)=sm2(n-1,(i-1)*1000+m-1)*exp(B(i-1)+as<double>(rnorm(1,0,sig[i-1])));
  }
  }
  }
  return wrap(sm2);'
  
  mcsim<-cxxfunction(signature(RB='vector',Rsig='vector',Rprice='vector'),
                     body=mcbody, plugin="Rcpp")
  mcs1=mcsim(B,sig,rep(1,50))
  kk=mcs1[61,]
  kkk=matrix(0,nrow=50,ncol=1000)
  for(i in 1:1000){
    cl=Css(i)
    kkk[,i]=kk[cl]
  }
  mcs1=kkk
  
  
  
  
  ##########################################################################################################
  
  ### risk sensitive portfolio optmizer. 
   # xi is risk tolerance parameter , IH=initial holding
  
  RSopt=function(xi,IH,price){
    
    eval_f <- function( x ) {RS_GBM(xi,x,price,mcs1)}
    
    eval_g_eq <- function( x ) {
      constr <- c( x%*%price -IH )
      grad <- c( price)
      return( list( "constraints"=constr, "jacobian"=grad ) )
    }
    
    # initial values and bounds
    lb <- -(0.1*IH)/abs(price)
    ub <- (0.1*IH)/abs(price)
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
  

  
  
  
###############################################################################################################
  
  #### calculate portfolio value of GBM model of risk sensitive reward function
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  RS_portfolio=function(xi,N){
    X=matrix(0,nrow=91,ncol=N)
    for(n in 1:N){
      
      ss=Css(n)#C++ 
      price=SimPrice[,ss]
      
      portfolio=matrix(0,nrow=3,ncol=20)
      portfolio[1,]=RSopt(xi,10000,In.price)$solution
      X[1:61,n]=price[1:61,]%*%portfolio[1,]
      for(i in 1:2){
        IH=X[30*i+1,n]
        c.price=price[30*i+1,]
        portfolio[i+1,]=RSopt(xi,IH,c.price)$solution
        X[(i*30+2):((i+1)*30+1),n]=price[(i*30+2):((i+1)*30+1),]%*%portfolio[i+1,]
      }
    }
    return(X)
  }
  
  
  ###############################################################################################################
  
  #### simulate price paths used in monte carlo method of optimization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  mmbody<-'
Environment stats("package:stats");
Function rnorm = stats["rnorm"];
NumericMatrix sm2(31,20000);

NumericVector B(RB);
NumericVector sig(Rsig);



for (int i=1; i<=1000; i++) {
for (int j=1; j<=20; j++) {
sm2(1-1,(i-1)*20+j-1)=1;
for (int m=1; m<=30; m++) {
sm2(m+1-1,(i-1)*20+j-1)=sm2(m-1,(i-1)*20+j-1)*exp(B(j-1)+as<double>(rnorm(1,0,sig[j-1])));
}
}
}
return wrap(sm2);'

Smm<-cxxfunction(signature(RB='vector',Rsig='vector'),
                 body=mmbody, plugin="Rcpp")
mxx<-Smm(B,sig)




###############################################################################################################

#### mean variance max min reward function and its gradient
## C++code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  maxmin<-'
double gamma  = as<double>(Rgamma);
double lambda = as<double>(Rlambda);
double theta = as<double>(Rtheta);

NumericVector x(Rx);
NumericVector E(RE);
NumericVector VA(RVA);
NumericVector K(RK);
NumericVector price(Rprice);
NumericMatrix maxs(Rmaxs);

NumericVector fun(1);
NumericVector funk(1);
NumericVector gk(1);
NumericMatrix Vs(31,20);


for(int i=1; i<=20; i++){
NumericVector p1(1);
for(int k=1; k<=20; k++){
p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1);
}
fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)+lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1));
}

NumericVector PX(1000);
for(int j=1; j<=1000; j++){
for(int k=1; k<=20; k++){
Vs(_,k-1)=maxs(_,(j-1)*20+k-1);
}
NumericVector XP(31);
for(int i=1; i<=31; i++){
for(int m=1; m<=20; m++){
XP(i-1)=XP(i-1)+x(m-1)*Vs(i-1,m-1)*price(m-1);}
}
PX(j-1)=max(XP)-min(XP);
}


funk=fun+theta*(sum(PX))*0.001;
return wrap(funk);
'
MaMi<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rtheta='numeric',
                               Rx='vector',RE='vector',RVA='vector',RK='vector',Rmaxs='Matrix',Rprice='Vector'),
                     body=maxmin, plugin="Rcpp")


###############################################################################################################

#optimizer of mean variance max min reward function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OP2=function(gamma,lambda,theta,IH,price){
  E=MCm2$E*price
  VA=MCm2$VA*(price^2)
  K=MCm2$K*(price^4)
  
  
  eval_f <- function(xx){
    MaMi(gamma,lambda,theta,xx,E,VA,K,mxx,price)
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    return( constr)
  }
  
  # initial values and bounds
  lb <- -(0.25*IH)/abs(price)
  ub <- (0.25*IH)/abs(price)
  
  x0 <- runif(20,lb,0.5/3*lb)
  x1 <- runif(20,0.5/3*lb,0.25/3*ub)
  x2 <- runif(20,0.25/3*ub,ub)
  opts = list(algorithm = "NLOPT_LN_COBYLA",
              "maxeval" = 100,
              xtol_rel = 1e-8)
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

########################################################
################################### portfolio value calculator

port.simM23=function(gamma,lambda,theta,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=20)
  portfolio[1,]=OP2(gamma,lambda,theta,1,In.price)$solution
  for(n in 1:N){
    ss=Css(n)#C++ 
    price=SimPrice[,ss]
    
    
    X[1:31,n]=price[1:31,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[30*i+1,n]
      c.price=price[30*i+1,]
      portfolio[i+1,]=OP2(gamma,lambda,theta,IH,c.price)$solution
      X[(i*30+2):((i+1)*30+1),n]=price[(i*30+2):((i+1)*30+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

  