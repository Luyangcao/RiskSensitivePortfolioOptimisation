library(PerformanceAnalytics)
library(nloptr)
library(Rcpp)
library(inline)
library(ggplot2)
ptm<-proc.time()
DDsim(rep(100,50),A,sd)-
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
  A=matrix(0,nrow=1,ncol=50)
  sd=matrix(0,nrow=1,ncol=50)
  
  
  for(i in 1:50){
    y=R.T[,i]
    x=R.T1[,i]
    L.Model=lm(y~x-1)
    A[1,i]=as.numeric(L.Model$coefficients)
    sd[1,i]=sd(L.Model$residuals)}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
}

########################################################################################################################

Sim_TruePath=function(){
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cxxs<-'NumericMatrix s(301,50000);
NumericVector A(RA);
NumericVector sd(Rsd);
NumericVector  price(Rprice);

for (int i=1; i<=50; i++) {
for(int m=1; m<=1000; m++){
s(1-1,(i-1)*1000+m-1)=price[i-1];
for(int n=1; n<=300; n++){
s(n+1-1,(i-1)*1000+m-1)=A[i-1]*s(n-1,(i-1)*1000+m-1)+as<double>(rnorm(1,0,sd[i-1]));
}
}
}
return wrap(s);'

Cpps<-cxxfunction(signature(RA='vector',Rsd='vector',Rprice='vector'),
                  body=cxxs, plugin="Rcpp")
s=Cpps(A,sd,In.price)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


}

########################################################################################################################

MC_sim=function(){

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mcsimbody<-'
Environment stats("package:stats");
Environment base("package:base");
Function mean = base["mean"];
Function var = stats["var"];
Function rnorm = stats["rnorm"];

NumericVector A(RA);
NumericVector sd(Rsd);

double IV = as<double>(RIV);
NumericMatrix pricesim(61,10000*50);



for(int i=1;i<=50;i++){
for(int m=1; m<=10000;m++){
pricesim(1-1,(i-1)*10000+m-1)=IV;
for(int n=1;n<=60;n++){
pricesim(n+1-1,(i-1)*10000+m-1)=A[i-1]*pricesim(n-1,(i-1)*10000+m-1)+as<double>(rnorm(1,0,sd[i-1]));
}
}}

NumericVector E(50);
NumericVector VA(50);
NumericVector K(50);


for(int i=1; i<=50;i++){
NumericVector st(10000);
for(int m=1; m<=10000;m++){
st[m-1]=pricesim(61-1,(i-1)*10000+m-1);}
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
ww<-cxxfunction(signature(RIV='numeric',RA='vector',Rsd='vector'),
                body=Mcsimbody, plugin="Rcpp")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MC<-ww(1,A,sd)


}

########################################################################################################################
%%%
EVALFF<-function(){
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
evalbody<-'
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
eval_fK<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rx='vector',RE='vector',RVA='vector',RK='vector'),
                     body=evalbody, plugin="Rcpp")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}


OPT2=function(gamma,lambda,IH,price){
  E=MC$E*price
  VA=MC$VA
  K=MC$K
  
  eval_f <- function( x ) {eval_fK(gamma,lambda,x,E,VA,K)}
  
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


port.sim2=function(gamma,lambda,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=s[,ss]
    
    portfolio=matrix(0,nrow=5,ncol=50)
    portfolio[1,]=OPT2(gamma,lambda,10000,In.price)$solution
    X[1:61,n]=price[1:61,]%*%portfolio[1,]
    for(i in 1:4){
      IH=X[60*i+1,n]
      c.price=price[60*i+1,]
      portfolio[i+1,]=OPT2(gamma,lambda,IH,c.price)$solution  
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

DDbody<-'
Environment stats("package:stats");
Environment base("package:base");
Function mean = base["mean"];
Function max = base["max"];
Function min = base["min"];
Function rnorm = stats["rnorm"];

NumericVector A(RA);
NumericVector sd(Rsd);
NumericVector Price(RPrice);


NumericMatrix pricesim(61,100*50);



for(int i=1;i<=50;i++){
  for(int m=1; m<=100;m++){
    pricesim(1-1,(i-1)*100+m-1)=Price(i-1);
    for(int n=1;n<=60;n++){
      pricesim(n+1-1,(i-1)*100+m-1)=A[i-1]*pricesim(n-1,(i-1)*100+m-1)+as<double>(rnorm(1,0,sd[i-1]));
    }
  }}

NumericVector DD(50);


for(int i=1; i<=50;i++){
  NumericVector Ndd(100);
   for(int j=1; j<=100;j++){
      NumericVector v(61);
      v=pricesim(_,((i-1)*100+j-1));
      Ndd(j-1)=as<double>(max(v))-as<double>(min(v));
   }
  DD(i-1)=as<double>(mean(Ndd));
}
 
return wrap(DD);
'
DDsim<-cxxfunction(signature(RPrice='vector',RA='vector',Rsd='vector'),
                body=DDbody, plugin="Rcpp")

MCDD1=DDsim(In.price,A,sd)
MCDD2=matrix(0,nrow=500,ncol=50)
MCDD3=matrix(0,nrow=500,ncol=50)
MCDD4=matrix(0,nrow=500,ncol=50)
MCDD5=matrix(0,nrow=500,ncol=50)
for(n in 1:1){
  ss=Css(n)
  price2=s[61,ss]
  price3=s[121,ss]
  price4=s[181,ss]
  price5=s[241,ss]
  MCDD2[n,]=DDsim(price2,A,sd)
  MCDD3[n,]=DDsim(price3,A,sd)
  MCDD4[n,]=DDsim(price4,A,sd)
  MCDD5[n,]=DDsim(price5,A,sd)
}

########################################################################################################################
OnlyGamma=function(){
  p.1000.0<-port.sim2(100,0,500)
  p.100.0<-port.sim2(10,0,500)
  p.20.0<-port.sim2(2,0,500)
  p.10.0<-port.sim2(1,0,500)
  p.5.0<-port.sim2(0.5,0,500)
  p.1.0<-port.sim2(0.1,0,500)
  p.02.0<-port.sim2(0.02,0,500)
  p.01.0<-port.sim2(0.01,0,500)
  p.005.0<-port.sim2(0.005,0,500)
  p.004.0<-port.sim2(0.004,0,500)
  p.003.0<-port.sim2(0.003,0,500)
  p.002.0<-port.sim2(0.002,0,500)
  p.001.0<-port.sim2(0.001,0,500)
  p.0005.0<-port.sim2(0.0005,0,500)
  p.0004.0<-port.sim2(0.0004,0,500)
  p.0001.0<-port.sim2(0.0001,0,500)
  p.00009.0<-port.sim2(9*0.1^5,0,500)
  p.00005.0<-port.sim2(5*0.1^5,0,500)
  p.00001.0<-port.sim2(0.1^5,0,500)
  
  

  f.1000.0<-perfYY(p.1000.0,500)
  f.100.0<-perfYY(p.100.0,500)
  f.20.0<-perfYY(p.20.0,500)
  f.10.0<-perfYY(p.10.0,500)
  f.5.0<-perfYY(p.5.0,500)
  f.1.0<-perfYY(p.1.0,500)
  f.02.0<-perfYY(p.02.0,500)
  f.01.0<-perfYY(p.01.0,500)
  f.005.0<-perfYY(p.005.0,500)
  f.004.0<-perfYY(p.004.0,500)
  f.003.0<-perfYY(p.003.0,500)
  f.002.0<-perfYY(p.002.0,500)
  f.001.0<-perfYY(p.001.0,500)
  f.0005.0<-perfYY(p.0005.0,500)
  f.0004.0<-perfYY(p.0004.0,500)
  f.0001.0<-perfYY(p.0001.0,500)
  f.00009.0<-perfYY(p.00009.0,500)
  f.00005.0<-perfYY(p.00005.0,500)
  f.00001.0<-perfYY(p.00001.0,500)
}
r00=rbind(f.1000.0,f.100.0,f.10.0,f.5.0,f.1.0,f.05.0,f.02.0,f.01.0,f.005.0,f.004.0,f.003.0,f.002.0,
          f.001.0,f.0005.0,f.0004.0,f.0001.0,f.00009.0,f.00005.0,f.00001.0)



Onlygammafigure=function(){
plot(density(p.100.0[301,]),xlim=c(10000,30000))
lines(density(p.10.0[301,]),col=2)
lines(density(p.1.0[301,]),col=4)
lines(density(p.02.0[301,]),col=5)
lines(density(p.01.0[301,]),col=5)
lines(density(p.005.0[301,]),col=6)
lines(density(p.004.0[301,]),col=7)
lines(density(p.003.0[301,]),col=8)
lines(density(p.002.0[301,]),col=9)
lines(density(p.001.0[301,]),col=10)
lines(density(p.0005.0[301,]),col=11)
lines(density(p.0001.0[301,]),col=12)
lines(density(p.00001.0[301,]),col=13)



x1<-c(p.100.0[301,],p.1.0[301,],p.02.0[301,],p.01.0[301,],p.005.0[301,],
      p.002.0[301,],p.001.0[301,],p.0005.0[301,],p.0001.0[301,],p.00001.0[301,])
x2=c(rep("10",500),rep("0.1",500),rep("0.02",500),rep("0.01",500),rep("0.005",500),
     rep("0.002",500),rep("0.001",500),rep("0.0005",500),rep("0.0001",500),rep("0.00001",500))
k<-data.frame(x1,x2)
colnames(k) <- c("Dollar", "Gamma")
p<-ggplot(data=k,aes(x=Dollar,color=Gamma))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
#p<-p+labs(title='Density of Final Wealth For Different Gamma')
library(plyr)
mu <- ddply(k, "Gamma", summarise, grp.mean=mean(Dollar))
head(mu)
p<-p+geom_vline(data=mu, aes(xintercept=grp.mean, color=Gamma),linetype="dashed")
p
}
########################################################################################################################

G0.002lambda=function(){
  p.002.100<-port.sim2(0.002,10,500)
  p.002.10<-port.sim2(0.002,1,500)
  p.002.1<-port.sim2(0.002,0.1,500)
  p.002.2<-port.sim2(0.002,0.1^2,500)
  p.002.3<-port.sim2(0.002,0.1^3,500)
  p.002.4<-port.sim2(0.002,0.1^4,500)
  p.002.5<-port.sim2(0.002,0.1^5,500)
  p.002.6<-port.sim2(0.002,0.1^6,500)
  p.002.7<-port.sim2(0.002,0.1^7,500)
  p.002.8<-port.sim2(0.002,0.1^8,500)
  p.002.9<-port.sim2(0.002,0.1^9,500)
  p.002.10<-port.sim2(0.002,0.1^10,500)
  p.002.11<-port.sim2(0.002,0.1^11,500)
  p.002.13<-port.sim2(0.002,0.1^13,500)
  p.002.15<-port.sim2(0.002,0.1^15,500)
  p.002.18<-port.sim2(0.002,0.1^18,500)
  p.002.23<-port.sim2(0.002,0.1^23,500)

  f.002.100<-perfYY(p.002.100,500)
  f.002.10<-perfYY(p.002.10,500)
  f.002.1<-perfYY(p.002.1,500)
  f.002.2<-perfYY(p.002.2,500)
  f.002.3<-perfYY(p.002.3,500)
  f.002.4<-perfYY(p.002.4,500)
  f.002.5<-perfYY(p.002.5,500)
  f.002.6<-perfYY(p.002.6,500)
  f.002.7<-perfYY(p.002.7,500)
  f.002.8<-perfYY(p.002.8,500)
  f.002.9<-perfYY(p.002.9,500)
  f.002.10<-perfYY(p.002.10,500)
  f.002.11<-perfYY(p.002.11,500)
  f.002.12<-perfYY(p.002.12,500)
  f.002.13<-perfYY(p.002.13,500)
  f.002.15<-perfYY(p.002.15,500)
  f.002.18<-perfYY(p.002.18,500)
  f.002.23<-perfYY(p.002.23,500)
}

r02=rbind(f.002.100,f.002.1,f.002.2,f.002.3,f.002.4,f.002.5,f.002.6,f.002.7,f.002.8,f.002.9,f.002.10,f.002.11
      ,f.002.13,f.002.15,f.002.18,f.002.23,f.002.0)

figure<-function(){
  x1<-c(p.100.0[301,],p.002.5[301,],p.002.6[301,],p.002.7[301,],p.002.8[301,],
        p.002.9[301,],p.002.10[301,],p.002.11[301,],p.002.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p3<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 18000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p3<-p3+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p3
}
########################################################################################################################

G0.5lambda=function(){
  p.5.100<-port.sim2(0.5,10,500)
  p.5.10<-port.sim2(0.5,1,500)
  p.5.1<-port.sim2(0.5,0.1,500)
  p.5.2<-port.sim2(0.5,0.1^2,500)
  p.5.3<-port.sim2(0.5,0.1^3,500)
  p.5.4<-port.sim2(0.5,0.1^4,500)
  p.5.5<-port.sim2(0.5,0.1^5,500)
  p.5.6<-port.sim2(0.5,0.1^6,500)
  p.5.7<-port.sim2(0.5,0.1^7,500)
  p.5.8<-port.sim2(0.5,0.1^8,500)
  p.5.9<-port.sim2(0.5,0.1^9,500)
  p.5.10<-port.sim2(0.5,0.1^10,500)
  p.5.11<-port.sim2(0.5,0.1^11,500)
  p.5.13<-port.sim2(0.5,0.1^13,500)
  p.5.15<-port.sim2(0.5,0.1^15,500)
  p.5.18<-port.sim2(0.5,0.1^18,500)
  p.5.23<-port.sim2(0.5,0.1^23,500)
  
  f.5.100<-perfYY(p.5.100,500)
  f.5.10<-perfYY(p.5.10,500)
  f.5.1<-perfYY(p.5.1,500)
  f.5.2<-perfYY(p.5.2,500)
  f.5.3<-perfYY(p.5.3,500)
  f.5.4<-perfYY(p.5.4,500)
  f.5.5<-perfYY(p.5.5,500)
  f.5.6<-perfYY(p.5.6,500)
  f.5.7<-perfYY(p.5.7,500)
  f.5.8<-perfYY(p.5.8,500)
  f.5.9<-perfYY(p.5.9,500)
  f.5.10<-perfYY(p.5.10,500)
  f.5.11<-perfYY(p.5.11,500)
  f.5.12<-perfYY(p.5.12,500)
  f.5.13<-perfYY(p.5.13,500)
  f.5.15<-perfYY(p.5.15,500)
  f.5.18<-perfYY(p.5.18,500)
  f.5.23<-perfYY(p.5.23,500)
}

r05=rbind(f.5.100,f.5.1,f.5.2,f.5.3,f.5.4,f.5.5,f.5.6,f.5.7,f.5.8,f.5.9,f.5.10,f.5.11
      ,f.5.13,f.5.15,f.5.18,f.5.23,f.5.0)

figure<-function(){
  x1<-c(p.100.0[301,],p.5.5[301,],p.5.6[301,],p.5.7[301,],p.5.8[301,],
        p.5.9[301,],p.5.10[301,],p.5.11[301,],p.5.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p1<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 14000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p1<-p1+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p1
}
########################################################################################################################

G0.0001lambda=function(){
  p.0001.100<-port.sim2(0.0001,10,500)
  p.0001.10<-port.sim2(0.0001,1,500)
  p.0001.1<-port.sim2(0.0001,0.1,500)
  p.0001.2<-port.sim2(0.0001,0.1^2,500)
  p.0001.3<-port.sim2(0.0001,0.1^3,500)
  p.0001.4<-port.sim2(0.0001,0.1^4,500)
  p.0001.5<-port.sim2(0.0001,0.1^5,500)
  p.0001.6<-port.sim2(0.0001,0.1^6,500)
  p.0001.7<-port.sim2(0.0001,0.1^7,500)
  p.0001.8<-port.sim2(0.0001,0.1^8,500)
  p.0001.9<-port.sim2(0.0001,0.1^9,500)
  p.0001.10<-port.sim2(0.0001,0.1^10,500)
  p.0001.11<-port.sim2(0.0001,0.1^11,500)
  p.0001.13<-port.sim2(0.0001,0.1^13,500)
  p.0001.15<-port.sim2(0.0001,0.1^15,500)
  p.0001.18<-port.sim2(0.0001,0.1^18,500)
  p.0001.23<-port.sim2(0.0001,0.1^23,500)
  
  f.0001.100<-perfYY(p.0001.100,500)
  f.0001.10<-perfYY(p.0001.10,500)
  f.0001.1<-perfYY(p.0001.1,500)
  f.0001.2<-perfYY(p.0001.2,500)
  f.0001.3<-perfYY(p.0001.3,500)
  f.0001.4<-perfYY(p.0001.4,500)
  f.0001.5<-perfYY(p.0001.5,500)
  f.0001.6<-perfYY(p.0001.6,500)
  f.0001.7<-perfYY(p.0001.7,500)
  f.0001.8<-perfYY(p.0001.8,500)
  f.0001.9<-perfYY(p.0001.9,500)
  f.0001.10<-perfYY(p.0001.10,500)
  f.0001.11<-perfYY(p.0001.11,500)
  f.0001.12<-perfYY(p.0001.12,500)
  f.0001.13<-perfYY(p.0001.13,500)
  f.0001.15<-perfYY(p.0001.15,500)
  f.0001.18<-perfYY(p.0001.18,500)
  f.0001.23<-perfYY(p.0001.23,500)
}
r01=rbind(f.0001.100,f.0001.1,f.0001.2,f.0001.3,f.0001.4,f.0001.5,f.0001.6,f.0001.7,
      f.0001.8,f.0001.9,f.0001.10,f.0001.11,f.0001.13,f.0001.15,f.0001.18,f.0001.23,f.0001.0)


figure0.0001<-function(){
  plot(density(p.100.0[301,]),xlim=c(10000,30000))
  lines(density(p.0001.5[301,]),col=2)
  lines(density(p.0001.6[301,]),col=3)
  lines(density(p.0001.7[301,]),col=4)
  lines(density(p.0001.8[301,]),col=3)
  lines(density(p.0001.9[301,]),col=3)
  lines(density(p.0001.10[301,]),col=5)
  lines(density(p.0001.11[301,]),col=6)
  lines(density(p.0001.13[301,]),col=4)
  lines(density(p.0001.0[301,]),col=2)
  
  x1<-c(p.100.0[301,],p.0001.5[301,],p.0001.6[301,],p.0001.7[301,],p.0001.8[301,],
        p.0001.9[301,],p.0001.10[301,],p.0001.11[301,],p.0001.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p2<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p2<-p2+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p2
}

########################################################################################################################

G0.0005lambda=function(){
  p.0005.100<-port.sim2(0.0005,10,500)
  p.0005.1<-port.sim2(0.0005,0.1,500)
  p.0005.4<-port.sim2(0.0005,0.1^4,500)
  p.0005.5<-port.sim2(0.0005,0.1^5,500)
  p.0005.6<-port.sim2(0.0005,0.1^6,500)
  p.0005.7<-port.sim2(0.0005,0.1^7,500)
  p.0005.85<-port.sim2(0.0005,5*0.1^8,500)
  p.0005.8<-port.sim2(0.0005,0.1^8,500)
  p.0005.95<-port.sim2(0.0005,5*0.1^9,500)
  p.0005.9<-port.sim2(0.0005,0.1^9,500)
  p.0005.10<-port.sim2(0.0005,0.1^10,500)
  p.0005.11<-port.sim2(0.0005,0.1^11,500)
  p.0005.13<-port.sim2(0.0005,0.1^13,500)
  p.0005.15<-port.sim2(0.0005,0.1^15,500)
  
  
  
  f.0005.100<-perfYY(p.0005.100,500)
  f.0005.1<-perfYY(p.0005.1,500)
  f.0005.4<-perfYY(p.0005.4,500)
  f.0005.5<-perfYY(p.0005.5,500)
  f.0005.6<-perfYY(p.0005.6,500)
  f.0005.7<-perfYY(p.0005.7,500)
  f.0005.85<-perfYY(p.0005.85,500)
  f.0005.8<-perfYY(p.0005.8,500)
  f.0005.95<-perfYY(p.0005.95,500)
  f.0005.9<-perfYY(p.0005.9,500)
  f.0005.10<-perfYY(p.0005.10,500)
  f.0005.11<-perfYY(p.0005.11,500)
  f.0005.13<-perfYY(p.0005.13,500)
  f.0005.15<-perfYY(p.0005.15,500)
}
r0005=rbind(f.0005.100,f.0005.1,f.0005.4,f.0005.5,f.0005.6,f.0005.7,f.0005.85,
            f.0005.8,f.0005.95,f.0005.9,f.0005.10,f.0005.11,f.0005.13,f.0005.15,f.0005.0)


figure0.0005<-function(){
  plot(density(p.100.0[301,]),xlim=c(10000,30000))
  lines(density(p.0001.5[301,]),col=2)
  lines(density(p.0001.6[301,]),col=3)
  lines(density(p.0001.7[301,]),col=4)
  lines(density(p.0001.8[301,]),col=3)
  lines(density(p.0001.9[301,]),col=3)
  lines(density(p.0001.10[301,]),col=5)
  lines(density(p.0001.11[301,]),col=6)
  lines(density(p.0001.13[301,]),col=4)
  lines(density(p.0001.0[301,]),col=2)
  
  x1<-c(p.100.0[301,],p.0005.5[301,],p.0005.6[301,],p.0005.7[301,],p.0005.8[301,],
        p.0005.9[301,],p.0005.10[301,],p.0005.11[301,],p.0005.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p0005<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p0005<-p0005+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p0005
}
