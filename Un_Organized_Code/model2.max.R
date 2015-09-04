mmbody<-'
Environment stats("package:stats");
Function rnorm = stats["rnorm"];
NumericMatrix sm2(61,50000);

NumericVector B(RB);
NumericVector sig(Rsig);
NumericVector price(Rprice);


for (int i=1; i<=1000; i++) {
for (int j=1; j<=50; j++) {
sm2(1-1,(i-1)*50+j-1)=price[j-1];
for (int m=1; m<=60; m++) {
sm2(m+1-1,(i-1)*50+j-1)=sm2(m-1,(i-1)*50+j-1)*B(j-1)+as<double>(rnorm(1,0,sig[j-1]));
}
}
}
return wrap(sm2);'

Smm<-cxxfunction(signature(RB='vector',Rsig='vector',Rprice='vector'),
                 body=mmbody, plugin="Rcpp")
mxx<-Smm(A,sd,In.price)




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
NumericMatrix Vs(61,50);


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
NumericVector XP(61);
for(int i=1; i<=61; i++){
for(int m=1; m<=50; m++){
XP(i-1)=XP(i-1)+x(m-1)*Vs(i-1,m-1);}
}
PX(j-1)=max(XP)-min(XP);
}


funk=fun+theta*(sum(PX))*0.001;
return wrap(funk);
'
MaMi<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rtheta='numeric',
                            Rx='vector',RE='vector',RVA='vector',RK='vector',Rmaxs='Matrix',Rprice='Vector'),
                  body=maxmin, plugin="Rcpp")



OP2=function(gamma,lambda,theta,IH,price){
  E=MC$E*price
  VA=MC$VA
  K=MC$K
  
  
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
              "maxeval" = 500,
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



port.simM23=function(gamma,lambda,theta,N){
  X=matrix(0,nrow=301,ncol=N)
  portfolio=matrix(0,nrow=5,ncol=50)
  portfolio[1,]=OP2(gamma,lambda,theta,10000,In.price)$solution
  for(n in 1:N){
    ss=Css(n)#C++ 
    price=s[,ss]
    
    
    X[1:61,n]=price[1:61,]%*%portfolio[1,]
    for(i in 1:4){
      IH=X[60*i+1,n]
      c.price=price[60*i+1,]
      portfolio[i+1,]=OP2(gamma,lambda,theta,IH,c.price)$solution
      X[(i*60+2):((i+1)*60+1),n]=price[(i*60+2):((i+1)*60+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}

port.simM23=function(gamma,lambda,theta,N){
  X=matrix(0,nrow=301,ncol=N)
  portfolio=OP2(gamma,lambda,theta,1,In.price)$solution

  for(n in 1:N){
    ss=Css(n)#C++ 
    price=s[,ss]
    
    
X[,n]=price%*%portfolio
  }
  return(X)
}


dl=OP2(0.0000001,0,0,1,In.price)
l
dl$solution%*%In.price

plot(density(port.simM23(0.00000001,0,0.1,1000)[301,]))


0.0000001*10000
m0<-perfYY(port.simM23(0.0000001,0,0,500),500)
m001<-perfYY(port.simM23(0.0000001,0,0.01,500),500)
m003<-perfYY(port.simM23(0.0000001,0,0.03,500),500)
m005<-perfYY(port.simM23(0.0000001,0,0.05,500),500)
m008<-perfYY(port.simM23(0.0000001,0,0.08,500),500)
m01<-perfYY(port.simM23(0.0000001,0,0.1,500),500)
m05<-perfYY(port.simM23(0.0000001,0,0.5,500),500)
m1<-perfYY(port.simM23(0.0000001,0,1,500),500)
m2<-perfYY(port.simM23(0.0000001,0,2,500),500)
m4<-perfYY(port.simM23(0.0000001,0,4,500),500)
m6<-perfYY(port.simM23(0.0000001,0,6,500),500)
m10<-perfYY(port.simM23(0.0000001,0,10,500),500)
Css
L=port.simM23(0.0000001,0,0,500)[301,]
mean(L)
sd(L)

rbind(m0,m05,m1,m2,m4,m6,m10)[,c(8,9,1,3,4,5,6)]

-0.001095055-0.0003749607
0.012428705
plot(density(rghyp(100000,best[[38]])-0.001470016))
lines(density(rghyp(100000,best[[46]])),col=2)


library(fGarch)
x.vec = as.vector(garchSim(garchSpec(rseed = 1985), n = N)[,1]
L=garchFit(~ arma(1,0)+garch(1,1), data = (rett), trace = FALSE)
K=garchSim(spec = S, n = 36, n.start = 100, extended = FALSE)
S=garchSpec(model = list(mu=-1.970e-03,ar = -1.867e-01, omega=8.220e-06,alpha = 5.755e-01, beta = 5.663e-01))
K
S=garchSpec(model = list(LIS))
LIS
LL=summary(L)
library(timeSeries)
LIS=coef(L)
plot(rett,type='l')
plot(K,type='l')  
