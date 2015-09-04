library(PerformanceAnalytics)
library(nloptr)
library(Rcpp)
library(inline)
library(ggplot2)
library(ghyp)
library(rugarch)



########################################################################################################################



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



########################################################################################################################

#  garch model fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=list() #4,9,15,18,28,29,31,46
BBB=BBlack[,-c(4,9,15,18,23,28,29,31,41,46)]
for(i in 21:40){
  y=BBB[1:100,i]
  rett=log(y+1)
  res[[i-20]]=ugarchfit(data =rett, spec = ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE)))
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




########################################################################################################################


#following code is used to simulate price from GRACH


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Simulator=function(){
  SIM=matrix(0,nrow=480,ncol=200000)
  for(i in 1:20){
    sim=ugarchsim(res[[i]], n.sim = 480, m.sim = 10000, n.start=100)
    JK=matrix(0,nrow=480,ncol=10000)
    for(j in 1:10000){
      JK[,j]=cumsum(fitted(sim)[,j])}
    SIM[,((i-1)*10000+1):(i*10000)]=exp(JK)
  }
  return(SIM)
}
Sim=Simulator()
Simprice=rbind(rep(1,20000),Sim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


}

########################################################################################################################



#Following code is used to find Mean Vairance and kurtosis for each stock at end of the period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

########################################################################################################################


#  a simple function calculate the number of columns, used to confirm column of each stock!

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


#### calculate portfolio value of GARCH model of mean-variance kurtoisis reward function

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

