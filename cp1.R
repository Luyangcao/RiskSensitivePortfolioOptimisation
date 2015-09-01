PRICE100 <- read.csv("~/Desktop/Project/prices/100PRICE.csv")
DATA.T1<-PRICE100[2817:3116,2:51] # Data t-1
DATA.T<-PRICE100[2818:3117,2:51] # Data t

R.T1<-PRICE100[2317:2816,2:51] # Data t-1
R.T<-PRICE100[2318:2817,2:51] # Data t
In.price=as.numeric(DATA.T1[1,])


#######################################################################################################################
##
## fit model X_t=AX_t-1+tau N(0,1)
##
#######################################################################################################################
A=matrix(0,nrow=1,ncol=50)
sd=matrix(0,nrow=1,ncol=50)


for(i in 1:50){
  y=R.T[,i]
  x=R.T1[,i]
  L.Model=lm(y~x-1)
  A[1,i]=as.numeric(L.Model$coefficients)
  sd[1,i]=sd(L.Model$residuals)}


s=matrix(0,nrow=301,ncol=2000*10) # 3min run time
for(i in 1:10){
  for(m in 1:2000){
    s[1,(i-1)*2000+m]=In.price[i] # initial price
    for(n in 1:300){
      s[n+1,(i-1)*2000+m]=A[i]*s[n,(i-1)*2000+m]+rnorm(1,0,sd[i])
    }
  }}

cxxfunction(signature())

src <- 'double x = as<double>(Rx);
        for (int i=-1; i<4; i++) 
             {x=x+i;}
        return wrap(x);' 


l <- cxxfunction(signature(Rx="numeric"),
                 body=src, plugin="Rcpp")

l(0)

library(Rcpp)
library(inline)
cxxs<-'NumericMatrix s(301,20000);
NumericVector A(RA);
NumericVector sd(Rsd);
NumericVector  price(Rprice);
for (int i=1; i<=10; i++) {
  for(int m=1; m<=2000; m++){
    s(1-1,(i-1)*2000+m-1)=price[i-1];
     for(int n=1; n<=300; n++){
        double x = as<double>(rnorm(1,A[i-1],1));
       s(n+1-1,(i-1)*2000+m-1)=A[i-1]*s(n-1,(i-1)*2000+m-1)+x;
    }
  }
}
return wrap(s);'

Cpps<-cxxfunction(signature(RA='vector',Rsd='vector',Rprice='vector'),
                  body=cxxs, plugin="Rcpp")

s=Cpps(A[1:10],sd[1:10],In.price[1:10])

s[22,1]

cxxs<-'NumericVector s(10);
NumericVector A(RA);
for (int i=0; i<10; i++){
double x = as<double>(rnorm(1,A[i],1));
s[i]=x;
}
return wrap(s);'

cxxs<-'NumericVector s(10);
NumericVector A(RA);
for (int i=0; i<10;i++){
s[i]=A[i];
}
return wrap(s);'


Cpps<-cxxfunction(signature(RA='vector'),
                  body=cxxs, plugin="Rcpp")
Cpps(A[1:10])
A[1:10]

?rnorm.h
A




J<-'
Environment stats("package:stats");
Function var = stats["var"];
Environment base("package:base");
Function mean = base["mean"];
double x = as<double>(Rx);
double y = as<double>(Ry);

NumericVector z(10);
NumericVector s(1);
z = rnorm(10,x,y);
s = l(0);
return wrap(s);'

Crnorm<-cxxfunction(signature(Rx='numeric',Ry='numeric'),
                  body=J, plugin="Rcpp")

?Rcpp
var(rnorm(10,2,3))
Crnorm(2,3)


?mean

ptm<-proc.time()
mean(k)
proc.time()-ptm
ptm<-proc.time()
Cvar(k)
proc.time()-ptm
Cpps(10,2,0.00005)




J<-'
NumericVector z(Rz);
NumericVector s(1);
NumericVector v(1);
for(int i=1; i<=100000000;i++){
s=s+z[i-1];
}
s=s/100000000;

return wrap(s);'

Cvar<-cxxfunction(signature(Rz='numeric'),
                    body=J, plugin="Rcpp")


Cvar(k)

?rnorm.h

MC.sim=function(IV){
  price.sim=matrix(0,nrow=101,ncol=20000*10)
  for(i in 1:10){
    for(m in 1:20000){
      price.sim[1,(i-1)*20000+m]=IV # initial price
      for(n in 1:100){
        price.sim[n+1,(i-1)*20000+m]=A[i]*price.sim[n,(i-1)*20000+m]+rnorm(1,0,sd[i])
      }
    }}
  
  E=matrix(0,nrow=1,ncol=10)
  VA=matrix(0,nrow=1,ncol=10)
  K=matrix(0,nrow=1,ncol=10)
  DD=matrix(0,nrow=1,ncol=10)
  
  for(i in 1:10){
    s.t=price.sim[101,((i-1)*20000+1):((i-1)*20000+20000)]
    E[i]=mean(s.t)
    VA[i]=var(s.t)
    K[i]=mean((s.t-E[i])^4)
    dd=c()
  }
  
  return(list('E'=E,'VA'=VA,'K'=K,'DD'=DD))}





Mcsimbody<-'
Environment stats("package:stats");
Environment base("package:base");
Function mean = base["mean"];
Function var = stats["var"];
Function rnorm = stats["rnorm"];

NumericVector A(RA);
NumericVector sd(Rsd);

double IV = as<double>(RIV);
NumericMatrix pricesim(31,10000*100);
for(int i=1;i<=100;i++){
  for(int m=1; m<=10000;m++){
   pricesim(1-1,(i-1)*10000+m-1)=IV;
    for(int n=1;n<=30;n++){
      pricesim(n+1-1,(i-1)*10000+m-1)=A[i-1]*pricesim(n-1,(i-1)*10000+m-1)+as<double>(rnorm(1,0,sd[i-1]));
  }
  }}

NumericVector E(100);
NumericVector VA(100);
NumericVector K(100);


for(int i=1; i<=100;i++){
NumericVector st(10000);
  for(int m=1; m<=10000;m++){
st[m-1]=pricesim(30,(i-1)*10000+m);}
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
MC<-ww(1,A,sd)
List ret; ret["x"] = xx; ret["y"] = yy;
return(ret);




#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


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


optbody='
Environment stats("package:stats");
Environment base("package:base");
Function function = base["function"];
Function var = stats["var"];
Function rnorm = stats["rnorm"];

NumericVector E(RE);
NumericVector VA(RVA);
NumericVector K(RK);






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

'

eval<-'
  NumericVector fun(1);
  for (int i=1; i<=100; i++){
    NumericVector p1(1);
    for(int k=1; k<=100,k++){
      p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1)
    }
    fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)
    +lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1))
  }
'


OPT1=function(gamma,IH,price){
  E=MC$E*price[1:4]
  VA=MC$VA[1:4]
  eval_f=function(x) eval_fK(gamma,x,E,VA)
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
OPT1(0.5,10000,In.price[1:4])
library(nloptr)


eval_f <- function( x,E,VA) {
  fun=0
  for(i in 1:4){
    fun=fun-x[i]*E[i]+gamma*x[i]^2*VA[i]
  }
  g=c()
  for(i in 1:4){
    g[i]=-E[i]+gamma*2*x[i]*VA[i]
  }
  return( list( "objective" =fun,
                "gradient" = g ) )
}


library(Rcpp)
library(inline)

evalbody<-'
double gamma = as<double>(Rgamma);
double lambda = as<double>(Rlambda);
NumericVector x(Rx);
NumericVector E(RE);
NumericVector VA(RVA);
NumericVector K(RK);


NumericVector fun(1);
NumericVector g(4);
NumericVector p2(1);
NumericVector gk(1);

for(int i=1; i<=4; i++){
 NumericVector p1(1);
  for(int k=1; k<=4; k++){
  p1=p1+x(i-1)*x(i-1)*VA(i-1)*x(k-1)*x(k-1)*VA(k-1);
  }
 fun=fun-x(i-1)*E(i-1)+gamma*x(i-1)*x(i-1)*VA(i-1)+lambda*x(i-1)*x(i-1)*x(i-1)*x(i-1)*K(i-1)+lambda*3*(p1-x(i-1)*x(i-1)*VA(i-1)*x(i-1)*x(i-1)*VA(i-1));
}

for(int i=1;i<=4;i++){
 p2=p2+x(i-1)*x(i-1)*VA(i-1);
}


for(int i=1; i<=4; i++){
  gk=2*x(i-1)*VA(i-1)*(p2-x(i-1)*x(i-1)*VA(i-1));
  g(i-1)=-E(i-1)+gamma*2*x(i-1)*VA(i-1)+lambda*4*x(i-1)*x(i-1)*x(i-1)*K(i-1)+6*lambda*gk(0);
}


List ret; ret["objective"] = fun; ret["gradient"] = g;
return wrap(ret);
'
eval_fK<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rx='vector',RE='vector',RVA='vector',RK='vector'),
                body=evalbody, plugin="Rcpp")

X=rep(1)
k=as.vector(array(rep(1,4),4))
eval_fK(0.5,0.5,k,k,k,k)
###
OPT1=function(gamma,lambda,IH,price){
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
                      "xtol_rel" = 1.0e-5 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-5,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)$solution
  return(res)}

ptm<-proc.time()
nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_eq=eval_g_eq,opts=opts)
proc.time()-ptm
ptm<-proc.time()
OPT2(0.5,0.5,10000,In.price)
proc.time()-ptm
gamma=0.5
lambda=0.5
wOPT2(0.5,0.5,10000,In.price[1:4])


port.sim2=function(gamma,lambda,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=s[,ss]
    
    portfolio=matrix(0,nrow=4,ncol=100)
    portfolio[1,]=OPT2(gamma,lambda,10000,In.price)$solution
    X[1:76,n]=price[1:76,]%*%portfolio[1,]
    for(i in 1:3){
      IH=X[75*i+1]
      c.price=price[75*i+1,]
      portfolio[i+1,]=OPT2(gamma,lambda,IH,c.price)$solution  
      X[(i*75+2):((i+1)*75+1),n]=price[(i*75+2):((i+1)*+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

ptm<-proc.time()
k=port.sim2(0.5,0.5,1)
proc.time()-ptm
k=port.sim2(0.5,0.5,1)
plot(k,type="l")

ptm<-proc.time()
for(i in 1:1000){
ss=c()
for(j in 1:100){
  ss[j]=(j-1)*1000+1
}}
proc.time()-ptm



ptm<-proc.time()
for(i in 1:1000){
Css(1)}
proc.time()-ptm


ssbody<-'
double n = as<double>(Rn);
NumericVector ss(100);
for(int i=1; i<=100; i++){
 ss(i-1)=(i-1)*1000+n;
}
return wrap(ss);
'
Css<-cxxfunction(signature(Rn='numeric'),
                     body=ssbody, plugin="Rcpp")
Css(1)



ssbody<-'
double n = as<double>(Rn);
NumericVector ss(100);
for(int i=1; i<=100; i++){
ss(i-1)=(i-1)*1000+n;
}
return wrap(ss);
'
Css<-cxxfunction(signature(Rn='numeric'),
                 body=ssbody, plugin="Rcpp")

Css(1)



kb<-'
Function ssf(CCs);
double n = as<double>(Rn);
NumericVector z(100);
z=ssf(1);
z=z(1)*n;
return wrap(z);
'
kk<-cxxfunction(signature(Rn='numeric',CCs='function'),
                 body=kb, plugin="Rcpp")

kk(3,Css)

dim(s[,Css(1)])


port.sim2=function(gamma,lambda,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=s[,ss]
    
    portfolio=matrix(0,nrow=4,ncol=100)
    portfolio[1,]=OPT2(gamma,lambda,10000,In.price)$solution
    X[1:76,n]=price[1:76,]%*%portfolio[1,]
    for(i in 1:3){
      IH=X[75*i+1]
      c.price=price[75*i+1,]
      portfolio[i+1,]=OPT2(gamma,lambda,IH,c.price)$solution  
      X[(i*75+2):((i+1)*75+1),n]=price[(i*75+2):((i+1)*+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}


Cport<-'
Environment base("package:base");
Function sum = base["sum"];

Function opt(Ropt);
Function Css(RCss);
double gamma = as<double>(Rgamma);
double lambda = as<double>(Rlambda);
double N = as<double>(RN);
NumericMatrix s(Rs);
NumericVector Inprice(RInprice);


NumericMatrix X(301,N);
NumericMatrix Pr(301,1);
NumericVector Portfolio(50);
NumericVector ss(50);

for(int i=1;i<=N;i++){
NumericMatrix Price(301,50);
ss=Css(i);
for(int k=1;k<=50;k++){
NumericMatrix::Column Pr=s(_,ss(k-1));
Price(_,k-1)=Pr;
}
Portfolio=opt(gamma,lambda,10000,Inprice);
for(int m=1;m<=61;m++){
    NumericVector p1(61);
    for(int n=1;n<=50;n++){
     p1(n-1)=Price(m-1,n-1)*Portfolio(n-1);
}
X(m-1,i-1)=as<double>(sum(p1));}

for(int j=1;j<=4;j++){
Portfolio=opt(gamma,lambda,X(60*j,i-1),Price(60*j,_));

for(int m=1;m<=60;m++){
    NumericVector p2(60);
for(int n=1;n<=50;n++){
p2(n-1)=Price(j*60+1+m-1,n-1)*Portfolio(n-1);
}
X(60*j+1+m-1,i-1)=as<double>(sum(p2));}
}
}
return wrap(X);
'



CPort<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',RN='numeric',RInprice='Vector',Rs='Matrix',Ropt='function',RCss='function'),
                body=Cport, plugin="Rcpp")
ptm<-proc.time()
CPort(0.5,0.5,10,In.price,s,OPT1,Css)
proc.time()-ptm
CPort(0.5,0.5,100,In.price,s,OPT1,Css)
ptm<-proc.time()
port.sim2(0.5,0.5,10)
proc.time()-ptm
port.sim2(0.5,0.5,100)



Cport<-'

Function opt(Ropt);
Function Css(RCss);
double gamma = as<double>(Rgamma);
double lambda = as<double>(Rlambda);
double N = as<double>(RN);
NumericMatrix s(Rs);
NumericVector Inprice(RInprice);


NumericMatrix X(301,N);
NumericMatrix Pr(301,1);
NumericVector Portfolio(50);
NumericVector ss(50);

for(int i=1;i<=N;i++){
NumericMatrix Price(301,50);
ss=Css(i);
for(int k=1;k<=50;k++){
NumericMatrix::Column Pr=s(_,ss(k-1));
Price(_,k-1)=Pr;
}
Portfolio=opt(gamma,lambda,10000,Inprice);
}
return wrap(X);
'

OPT1(0.5,0.5,10000,In.price)
CPort<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',RN='numeric',RInprice='Vector',Rs='Matrix',Ropt='function',RCss='function'),
                   body=Cport, plugin="Rcpp")

CPort(0.5,0.5,1,In.price,s,OPT1,Css)


