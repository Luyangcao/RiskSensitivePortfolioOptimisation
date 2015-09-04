PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
DATA.T1<-PRICE100[3027:3116,21:40] # Data t-1
DATA.T<-PRICE100[3028:3117,21:40] # Data t

R.T1<-PRICE100[2317:3026,21:40] # Data t-1
R.T<-PRICE100[2318:3027,21:40] # Data t
In.price=as.numeric(DATA.T1[1,])

dbest=list()


for(i in 1:20){
  y=R.T[,i]
  x=R.T1[,i]
  rett<-log(y/x)
  fit <- stepAIC.ghyp(rett, dist = c("hyp", "NIG", "t","gauss"),symmetric = NULL, control = list(maxit = 5000),na.rm = T,silent = TRUE)
  bestfit <-fit$best.model
  dbest[i]=bestfit}


mtcxxsm2<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
NumericMatrix sm2(91,300000);

List best(Rbest);
NumericVector  price(Rprice);

for (int i=1; i<=20; i++) {
for(int m=1; m<=10000; m++){
sm2(1-1,(i-1)*10000+m-1)=price[i-1];
NumericVector  gg(90);
gg=rghyp(90,best[i-1]);
for(int n=1; n<=90; n++){
sm2(n+1-1,(i-1)*10000+m-1)=sm2(n-1,(i-1)*10000+m-1)*exp(gg[n-1]);
}
}
}
return wrap(sm2);'

mtCppsm2<-cxxfunction(signature(Rbest='List',Rprice='vector'),
                     body=mtcxxsm2, plugin="Rcpp")


ptm<-proc.time()
mtsm2W=mtCppsm2(dbest,In.price)
proc.time()-ptm



mtMcsimbodym2<-'
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
NumericMatrix pricesim(31,10000*20);



for(int i=1;i<=20;i++){
for(int m=1; m<=10000;m++){
pricesim(1-1,(i-1)*10000+m-1)=IV;
NumericVector  gg(30);
gg=rghyp(30,best[i-1]);
for(int n=1;n<=30;n++){
pricesim(n+1-1,(i-1)*10000+m-1)=pricesim(n-1,(i-1)*10000+m-1)*exp(gg[n-1]);
}
}}

NumericVector E(20);
NumericVector VA(20);
NumericVector K(20);
NumericVector v(30);

for(int i=1; i<=20;i++){
NumericVector st(10000);
for(int m=1; m<=10000;m++){
st[m-1]=pricesim(11-1,(i-1)*10000+m-1);
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
mtwwm2<-cxxfunction(signature(RIV='numeric',Rbest='List'),
                   body=mtMcsimbodym2, plugin="Rcpp")
mtMCm2<-mtwwm2(1,dbest)
mtMCm21<-mtwwm2(1,dbest)
mtMCm22<-mtwwm2(1,dbest)
mtMCm2.100<-mtwwm2(100,dbest)




mevalbody4<-'
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
meval_f4<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rx='vector',RE='vector',RVA='vector',RK='vector'),
                     body=mevalbody4, plugin="Rcpp")




OPTM3=function(gamma,lambda,IH,price){
  
  E=mtMCm2$E*price
  VA=mtMCm2$VA*price^2
  K=mtMCm2$K*price^4
  
  eval_f <- function( x ) {meval_f4(gamma,lambda,x,E,VA,K)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.25*IH)/abs(price)
  ub <- (0.25*IH)/abs(price)
  set.seed(998)
  x0 <- runif(20,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-8 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-8,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}


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



port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=20)
  portfolio[1,]=OPTM3(gamma,lambda,1,In.price)$solution
  weight=portfolio[1,]*In.price/1
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=mtsm2W[,ss]
    
    X[1:31,n]=price[1:31,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[30*i+1,n]
      c.price=price[30*i+1,]
      portfolio[i+1,]=weight*IH/c.price
      X[(i*30+2):((i+1)*30+1),n]=price[(i*30+2):((i+1)*30+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=OPTM3(gamma,lambda,1,In.price)$solution
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=mtsm2W[,ss]
    
    X[,n]=price%*%portfolio
  }
  return(X)
}



tevalbodyf5<-'
#include <math.h>

double xi = as<double>(Rxi);


NumericVector x(Rx);
NumericVector price(Rprice);
NumericMatrix mcs(Rmcs1);

NumericVector funk(1000);
NumericVector fun(1);
NumericVector g(20);
NumericVector Vmcs(20);
NumericVector XP(1000);

for(int j=1; j<=1000; j++){
Vmcs=mcs(_,j-1);
for(int i=1; i<=20; i++){
XP(j-1)=XP(j-1)+x(i-1)*Vmcs(i-1)*price(i-1);
}
funk(j-1)=exp(-xi*XP(j-1));
}
fun=log(sum(funk));


for(int k=1; k<=20; k++){
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



mtmcbody<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
NumericMatrix sm2(31,200000);

List best(Rbest);
NumericVector  price(Rprice);

for (int i=1; i<=20; i++) {
for(int m=1; m<=10000; m++){
sm2(1-1,(i-1)*10000+m-1)=price[i-1];
NumericVector  gg(30);
gg=rghyp(30,best[i-1]);
for(int n=1; n<=30; n++){
sm2(n+1-1,(i-1)*10000+m-1)=sm2(n-1,(i-1)*10000+m-1)*exp(gg[n-1]);
}
}
}
return wrap(sm2);'

mtmcsim<-cxxfunction(signature(Rbest='List',Rprice='vector'),
                    body=mtmcbody, plugin="Rcpp")

mtmcs1=mtmcsim(dbest,rep(1,20))
kk=mtmcs1[31,]
kkk=matrix(0,nrow=20,ncol=10000)
for(i in 1:10000){
  cl=Css(i)
  kkk[,i]=kk[cl]
}
mtmcs1=kkk



OPTM25=function(xi,IH,price){
  
  eval_f <- function( x ) {teval_f5(xi,x,price,mtmcs1)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.25*IH)/abs(price)
  ub <- (0.25*IH)/abs(price)
  set.seed(998)
  x0 <- runif(20,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-8 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-8,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
  return(res)}




port.simM25=function(xi,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=20)
  portfolio[1,]=OPTM25(xi,1,In.price)$solution
  weight=portfolio[1,]*In.price/1
  for(n in 1:N){
    ss=Css(n)#C++ 
    price=mtsm2W[,ss]
    
    
    X[1:31,n]=price[1:31,]%*%portfolio[1,]
    for(i in 1:2){
      IH=X[30*i+1,n]
      c.price=price[30*i+1,]
      portfolio[i+1,]=weight*IH/c.price
      X[(i*30+2):((i+1)*30+1),n]=price[(i*30+2):((i+1)*30+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}


port.simM25=function(xi,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=OPTM25(xi,1,In.price)$solution
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=mtsm2W[,ss]
    
    X[,n]=price%*%portfolio
  }
  return(X)
}


perfYY=function(YYY,N){
  X=YYY
  SR=numeric(N);MDD=numeric(N);
  for(i in 1:N){
    k=X[,i]
    ret=diff(k)/k[1:90] #return
    R=data.frame(ret,row.names=as.Date(z))
    SR[i]=SharpeRatio(R,Rb,FUN="StdDev") # 
    MDD[i]=maxDrawdown(R)#
  }
  return(list('SR'=mean(SR),'MDD'=mean(MDD)))}