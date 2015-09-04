library(PerformanceAnalytics)
library(nloptr)
library(Rcpp)
library(inline)
library(ggplot2)
library(ghyp)
library(rugarch)



########################################################################################################################



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



########################################################################################################################

## fit GH distribution best fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best=list()

for(i in 1:50){
  y=BBlack[1:72,1]
  rett<-log(y+1)
  fit <- stepAIC.ghyp(rett, dist = c( "hyp", "NIG", "t", "gauss"),symmetric = NULL, control = list(maxit = 5000),na.rm = T,silent = TRUE)
  bestfit <-fit$best.model
  best[i]=bestfit}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




########################################################################################################################


#following code is used to simulate price from GH model.C++ code


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PriceSimbody<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
NumericMatrix sm2(37,500000);

List best(Rbest);
NumericVector  price(Rprice);

for (int i=1; i<=50; i++) {
for(int m=1; m<=10000; m++){
sm2(1-1,(i-1)*10000+m-1)=price[i-1];
NumericVector  gg(36);
gg=rghyp(36,best[i-1]);
for(int n=1; n<=36; n++){
sm2(n+1-1,(i-1)*10000+m-1)=sm2(n-1,(i-1)*10000+m-1)*exp(gg[n-1]);
}
}
}
return wrap(sm2);'

PriceSim<-cxxfunction(signature(Rbest='List',Rprice='vector'),
                      body=tcxxsm2, plugin="Rcpp")

SimPrice=PriceSim(best,rep(1,50))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


}

########################################################################################################################



#Following code is used to find Mean Vairance and kurtosis for each stock at end of the period. C++ code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MCsimbody<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
Environment stats("package:stats");
Environment base("package:base");

Function mean = base["mean"];
Function var = stats["var"];
Function rt = stats["rt"];
Function max = base["max"];
Function min = base["min"];

List best(Rbest);

double IV = as<double>(RIV);
NumericMatrix pricesim(13,100000*50);



for(int i=1;i<=50;i++){
for(int m=1; m<=100000;m++){
pricesim(1-1,(i-1)*100000+m-1)=IV;
NumericVector  gg(12);
gg=rghyp(12,best[i-1]);
for(int n=1;n<=12;n++){
pricesim(n+1-1,(i-1)*100000+m-1)=pricesim(n-1,(i-1)*100000+m-1)*exp(gg[n-1]);
}
}}

NumericVector E(50);
NumericVector VA(50);
NumericVector K(50);


for(int i=1; i<=50;i++){
NumericVector st(100000);
for(int m=1; m<=100000;m++){
st[m-1]=pricesim(13-1,(i-1)*100000+m-1);
}
E[i-1]=as<double>(mean(st));
VA[i-1]=as<double>(var(st));
NumericVector st2(100000);
for(int n=1; n<=100000;n++){
st2[n-1]=(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1]);}
K[i-1]=as<double>(mean(st2));
}
List ret; ret["E"] = E; ret["VA"] = VA;ret["K"] = K;
return wrap(ret);
'
MCproperty<-cxxfunction(signature(RIV='numeric',Rbest='List'),
                        body=MCsimbody, plugin="Rcpp")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MCm2<-MCproperty(1,B,sig)
tMCm2<-twwm2(1,best)



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
NumericVector g(50);
NumericVector p2(1);
NumericVector gk(1);

for(int i=1; i<=50; i++){
NumericVector p1(1);
for(int k=1; k<=50; k++){
p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1);
}
fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)+lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1));
}

for(int i=1;i<=50;i++){
p2=p2+x(i-1)*x(i-1)*VA(i-1);
}


for(int i=1; i<=50; i++){
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

OPTM3=function(gamma,lambda,IH,price){    
  E=tMCm2$E*price
  VA=tMCm2$VA*(price^2)
  K=tMCm2$K*(price^4)
  
  eval_f <- function( x ) {eval_f4(gamma,lambda,x,E,VA,K)}
  
  eval_g_eq <- function( x ) {    # equality constriant
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.1*IH)/abs(price)
  ub <- (0.1*IH)/abs(price)
  set.seed(998)
  x0 <- runif(50,lb,0.5/3*lb)
  x1 <- runif(50,0.5/3*lb,0.25/3*ub)
  x2 <- runif(50,0.25/3*ub,ub)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssbody<-'
double n = as<double>(Rn);
NumericVector ss(50);
for(int i=1; i<=50; i++){
ss(i-1)=(i-1)*10000+n;
}
return wrap(ss);
'
Css<-cxxfunction(signature(Rn='numeric'),
                 body=ssbody, plugin="Rcpp")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#### calculate portfolio value of GBM model of mean-variance kurtoisis reward function
port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=37,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=50)
  portfolio[1,]=OPTM3(gamma,lambda,1,rep(1,50))$solution
  weight=portfolio[1,]*rep(1,50)/1
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=tsm2W[,ss]
    
    X[1:13,n]=price[1:13,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[12*i+1,n]
      c.price=price[12*i+1,]
      portfolio[i+1,]=weight*IH/c.price
      X[(i*12+2):((i+1)*12+1),n]=price[(i*12+2):((i+1)*12+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

########################################################################################################################



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RF=rep(0,50)
z <- as.yearmon(1960 + seq(0, 11)/12)
Rf=data.frame(RF,row.names=as.Date(z))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# calcualte the perfomance indices, SR=sharpe ratio, IR=informantion ratio, ADD=average drawdown, MDD=maximum drawdown
# NDD=number of drawdowan, MT=maximum time of drawdown, LS=loss probability, Re=return, V=volatility

perfYY=function(YYY,N){
  X=YYY
  SR=numeric(N);IR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N);MT=numeric(N);LS=numeric(N);Re=numeric(N)
  for(i in 1:N){
    k=X[,i]
    ret=diff(k)/k[1:36] #return
    R=data.frame(ret,row.names=as.Date(z))
    SR[i]=SharpeRatio(R,Rf,FUN="StdDev") # 
    ADD[i]=AverageDrawdown(R) #
    MDD[i]=maxDrawdown(R)#
    FD=findDrawdowns(R)
    NDD[i]=length(FD$from)
    MT[i]=max(FD$to-FD$from)
    LS[i]=as.numeric(k[37]<1)
    Re[i]=k[37]-1
    IR[i]=InformationRatio(R,Rf)
  }
  return(list('SR'=mean(SR),'IR'=mean(IR),'MDD'=mean(MDD),'ADD'=mean(ADD),'NDD'=mean(NDD),
              "MT"=mean(MT),"LS"=mean(LS),"R"=mean(Re),"V"=sd(Re)))}


########################################################################################################################

########################################################################################################################



###risk sensitive reward function and its gradient, C++ code

## mcs or mcs1 is simulated terminal value of each stock, used to monte carlo method in optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tevalbodyf5<-'
#include <math.h>

double xi = as<double>(Rxi);


NumericVector x(Rx);
NumericVector price(Rprice);
NumericMatrix mcs(Rmcs1);

NumericVector funk(1000);
NumericVector fun(1);
NumericVector g(50);
NumericVector Vmcs(50);
NumericVector XP(1000);

for(int j=1; j<=1000; j++){
Vmcs=mcs(_,j-1);
for(int i=1; i<=50; i++){
XP(j-1)=XP(j-1)+x(i-1)*Vmcs(i-1)*price(i-1);
}
funk(j-1)=exp(-xi*XP(j-1));
}
fun=log(sum(funk));


for(int k=1; k<=50; k++){
NumericVector gfz(1000);
for(int j=1; j<=1000; j++){
Vmcs=mcs(_,j-1);
gfz[j-1]=exp(-xi*XP(j-1))*(-xi*Vmcs[k-1]*price(k-1));}
g(k-1)=sum(gfz)/(sum(funk));
}


List ret; ret["objective"] = fun; ret["gradient"] = g;
return wrap(ret);
'
teval_f5<-cxxfunction(signature(Rxi='numeric',Rx='vector',Rprice='vector',Rmcs1='Matrix'),
                      body=tevalbodyf5, plugin="Rcpp")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






## tmsim is simulated terminal value of each stock, used to monte carlo method in optimization

tmcbody<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
NumericMatrix sm2(13,500000);

List best(Rbest);
NumericVector  price(Rprice);

for (int i=1; i<=50; i++) {
for(int m=1; m<=10000; m++){
sm2(1-1,(i-1)*10000+m-1)=price[i-1];
NumericVector  gg(12);
gg=rghyp(12,best[i-1]);
for(int n=1; n<=12; n++){
sm2(n+1-1,(i-1)*10000+m-1)=sm2(n-1,(i-1)*10000+m-1)*exp(gg[n-1]);
}
}
}
return wrap(sm2);'

tmcsim<-cxxfunction(signature(Rbest='List',Rprice='vector'),
                    body=tmcbody, plugin="Rcpp")

tmcs1=tmcsim(best,rep(1,50))
kk=tmcs1[13,]
kkk=matrix(0,nrow=50,ncol=1000)
for(i in 1:1000){
  cl=Css(i)
  kkk[,i]=kk[cl]
}
tmcs1=kkk



##########################################################################################################

### risk sensitive portfolio optmizer. 
# xi is risk tolerance parameter , IH=initial holding

RSopt=function(xi,IH,price){
  
  eval_f <- function( x ) {tevalbodyf5(xi,x,price,mcs1)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.1*IH)/abs(price)
  ub <- (0.1*IH)/abs(price)
  set.seed(998)
  x0 <- runif(50,lb,0.5/3*lb)
  x1 <- runif(50,0.5/3*lb,0.25/3*ub)
  x2 <- runif(50,0.25/3*ub,ub)
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
  
  RS_GH_port=function(xi,N){
    X=matrix(0,nrow=37,ncol=N)
    portfolio=matrix(0,nrow=3,ncol=50)
    portfolio[1,]=OPTM25(xi,1,rep(1,50))$solution
    weight=portfolio[1,]*rep(1,50)/1
    for(n in 1:N){
      ss=Css(n)#C++ 
      price=SimPrice[,ss]
      
      
      X[1:13,n]=price[1:13,]%*%portfolio[1,]
      for(i in 1:2){
        IH=X[12*i+1,n]
        c.price=price[12*i+1,]
        portfolio[i+1,]=weight*IH/c.price
        X[(i*12+2):((i+1)*12+1),n]=price[(i*12+2):((i+1)*12+1),]%*%portfolio[i+1,]
      }
    }
    return(X)
  }


###############################################################################################################

#### simulate price paths used in monte carlo method of optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  mmbody<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
NumericMatrix sm2(13,50000);

NumericVector best(R);




for (int i=1; i<=1000; i++) {
for (int j=1; j<=12; j++) {
sm2(1-1,(i-1)*50+j-1)=1;
for (int m=1; m<=12; m++) {
sm2(m+1-1,(i-1)*50+j-1)=sm2(m-1,(i-1)*50+j-1)*exp(as<double>(rghyp(1,best[i-1])));
}
}
}
return wrap(sm2);'

Smm<-cxxfunction(signature(Rbest='list'),
                 body=mmbody, plugin="Rcpp")
mxx<-Smm(best)




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
NumericMatrix Vs(13,50);


for(int i=1; i<=50; i++){
NumericVector p1(1);
for(int k=1; k<=50; k++){
p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1);
}
fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)+lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1));
}

NumericVector PX(1000);
for(int j=1; j<=1000; j++){
for(int k=1; k<=50; k++){
Vs(_,k-1)=maxs(_,(j-1)*50+k-1);
}
NumericVector XP(13);
for(int i=1; i<=13; i++){
for(int m=1; m<=50; m++){
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
  E=tMCm2$E*price
  VA=tMCm2$VA*(price^2)
  K=tMCm2$K*(price^4)
  
  
  eval_f <- function(xx){
    MaMi(gamma,lambda,theta,xx,E,VA,K,mxx,price)
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    return( constr)
  }
  
  # initial values and bounds
  lb <- -(0.1*IH)/abs(price)
  ub <- (0.1*IH)/abs(price)
  
  x0 <- runif(50,lb,0.5/3*lb)
  x1 <- runif(50,0.5/3*lb,0.25/3*ub)
  x2 <- runif(50,0.25/3*ub,ub)
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
  X=matrix(0,nrow=37,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=50)
  portfolio[1,]=OP2(gamma,lambda,theta,1,In.price)$solution
  for(n in 1:N){
    ss=Css(n)#C++ 
    price=SimPrice[,ss]
    
    
    X[1:13,n]=price[1:13,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[12*i+1,n]
      c.price=price[12*i+1,]
      portfolio[i+1,]=OP2(gamma,lambda,theta,IH,c.price)$solution
      X[(i*12+2):((i+1)*12+1),n]=price[(i*12+2):((i+1)*12+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

