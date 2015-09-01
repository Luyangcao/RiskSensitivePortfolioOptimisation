PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
DATA.T1<-PRICE100[2817:3116,2:101] # Data t-1
DATA.T<-PRICE100[2818:3117,2:101] # Data t

R.T1<-PRICE100[2317:2816,2:101] # Data t-1
R.T<-PRICE100[2318:2817,2:101] # Data t
In.price=as.numeric(DATA.T1[1,])

ptm<-proc.time()
OPT2(0.5,0.1,10000,In.price)
proc.time()-ptm
#######################################################################################################################
##
## fit model X_t=AX_t-1+tau N(0,1)
##
#######################################################################################################################
A=matrix(0,nrow=1,ncol=100)
sd=matrix(0,nrow=1,ncol=100)
library(Rcpp)
library(inline)
library(nloptr)
for(i in 1:100){
  y=R.T[,i]
  x=R.T1[,i]
  L.Model=lm(y~x-1)
  A[1,i]=as.numeric(L.Model$coefficients)
  sd[1,i]=sd(L.Model$residuals)}



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

ptm<-proc.time()
s=Cpps(A,sd,In.price)
proc.time()-ptm


########M  MC.sim
for(){
  ############################################################
  ############################################################
  ############################################################
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
##########################################################################################
##########################################################################################
##########################################################################################
ptm<-proc.time()
MC<-ww(1,A,sd)
proc.time()-ptm


}




#########eval_f
for(){
  ######################################################################## 
  ######################################################################## 
  ######################################################################## 
  ######################################################################## 
  
  
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

  ######################################################################## 
  ######################################################################## 
  ######################################################################## 
  ######################################################################## 
  ######################################################################## 
  
  
  
  
}

############
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

OPT2(0.5,0.1,10000,In.price)
########## Css function

for(){
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

}
#########################################################################################################################

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


plot(density(p.20.0[301,]))

#########################################################################################################################

ptm<-proc.time()
p.20.0<-port.sim2(2,0,300)
proc.time()-ptm
000=function(){
  p.1000.0<-port.sim2(100,0,300)
  p.200.0<-port.sim2(20,0,300)
  p.100.0<-port.sim2(10,0,300)
  p.20.0<-port.sim2(2,0,300)
  p.10.0<-port.sim2(1,0,300)
  p.5.0<-port.sim2(0.5,0,300)
  p.1.0<-port.sim2(0.1,0,300)
  p.01.0<-port.sim2(0.01,0,300)
  p.005.0<-port.sim2(0.005,0,300)
  p.004.0<-port.sim2(0.004,0,300)
  p.003.0<-port.sim2(0.003,0,300)
  p.002.0<-port.sim2(0.002,0,300)
  p.001.0<-port.sim2(0.001,0,300)
  p.0008.0<-port.sim2(0.0008,0,300)
  p.0007.0<-port.sim2(0.0007,0,300)
  p.0006.0<-port.sim2(0.0006,0,300)
  p.0005.0<-port.sim2(0.0005,0,300)
  p.0004.0<-port.sim2(0.0004,0,300)
  p.0001.0<-port.sim2(0.0001,0,300)
  
  f.1000.0<-perfYY(p.1000.0,300)
  f.200.0<-perfYY(p.200.0,300)
  f.100.0<-perfYY(p.100.0,300)
  f.20.0<-perfYY(p.20.0,300)
  f.10.0<-perfYY(p.10.0,300)
  f.5.0<-perfYY(p.5.0,300)
  f.1.0<-perfYY(p.1.0,300)
  f.01.0<-perfYY(p.01.0,300)
  f.005.0<-perfYY(p.005.0,300)
  f.004.0<-perfYY(p.004.0,300)
  f.003.0<-perfYY(p.003.0,300)
  f.002.0<-perfYY(p.002.0,300)
  f.001.0<-perfYY(p.001.0,300)
  f.0008.0<-perfYY(p.0008.0,300)
  f.0007.0<-perfYY(p.0007.0,300)
  f.0006.0<-perfYY(p.0006.0,300)
  f.0005.0<-perfYY(p.0005.0,300)
  f.0004.0<-perfYY(p.0004.0,300)
  f.0001.0<-perfYY(p.0001.0,300)
}


rbind(f.1000.0,f.200.0,f.100.0,f.20.0,f.10.0,f.5.0,f.1.0,f.01.0,f.005.0,f.004.0,f.003.0,f.002.0,
      f.001.0,f.0008.0,f.0007.0,f.0006.0,
      f.0005.0,f.0004.0,f.0001.0)


plot(density(p.20.0[301,]),xlim=c(10000,30000))
lines(density(p.100.0[301,]),col=2)
lines(density(p.200.0[301,]),col=3)
lines(density(p.1000.0[301,]),col=4)
lines(density(p.10.0[301,]),col=2)
lines(density(p.5.0[301,]),col=3)
lines(density(p.1.0[301,]),col=4)
lines(density(p.01.0[301,]),col=4)
lines(density(p.005.0[301,]),col=4)
lines(density(p.003.0[301,]),col=4)
lines(density(p.001.0[301,]),col=4)
lines(density(p.0008.0[301,]),col=4)
lines(density(p.0005.0[301,]),col=4)
lines(density(p.0001.0[301,]),col=4)
#########################################################################################################################
G0.01=function(){
  p.01.20<-port.sim2(0.01,0.20,300)
  p.01.1<-port.sim2(0.01,0.1,300)
  p.01.2<-port.sim2(0.01,0.1^2,300)
  p.01.3<-port.sim2(0.01,0.1^3,300)
  p.01.4<-port.sim2(0.01,0.1^4,300)
  p.01.5<-port.sim2(0.01,0.1^5,300)
  p.01.6<-port.sim2(0.01,0.1^6,300)
  p.01.7<-port.sim2(0.01,0.1^7,300)
  p.01.8<-port.sim2(0.01,0.1^8,300)#
  p.01.9<-port.sim2(0.01,0.1^9,300)
  p.01.10<-port.sim2(0.01,0.1^10,300)
  p.01.11<-port.sim2(0.01,0.1^11,300)
  p.01.12<-port.sim2(0.01,0.1^12,300)
  p.01.13<-port.sim2(0.01,0.1^13,300)
  p.01.14<-port.sim2(0.01,0.1^14,300)
  p.01.15<-port.sim2(0.01,0.1^15,300)
  p.01.18<-port.sim2(0.01,0.1^18,300)
  p.01.22<-port.sim2(0.01,0.1^22,300)
  
  
  
  
  f.01.20<-perfYY(p.01.20,300)
  f.01.1<-perfYY(p.01.1,300)
  f.01.2<-perfYY(p.01.2,300)
  f.01.3<-perfYY(p.01.3,300)
  f.01.4<-perfYY(p.01.4,300)
  f.01.5<-perfYY(p.01.5,300)
  f.01.6<-perfYY(p.01.6,300)
  f.01.7<-perfYY(p.01.7,300)
  f.01.8<-perfYY(p.01.8,300)
  f.01.9<-perfYY(p.01.9,300)
  f.01.10<-perfYY(p.01.10,300)
  f.01.11<-perfYY(p.01.11,300)
  f.01.12<-perfYY(p.01.12,300)
  f.01.13<-perfYY(p.01.13,300)
  f.01.14<-perfYY(p.01.14,300)
  f.01.15<-perfYY(p.01.15,300)
  f.01.18<-perfYY(p.01.18,300)
  f.01.22<-perfYY(p.01.22,300)
}

rbind(f.01.20,f.01.1,f.01.2,f.01.3,f.01.4,f.01.5,f.01.6,f.01.7,f.01.8,f.01.9,f.01.10,
      f.01.11,f.01.12,f.01.13,f.01.14,f.01.15,f.01.18,f.01.0)##
