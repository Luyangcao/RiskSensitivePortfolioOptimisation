

ptm<-proc.time()
port.simM25(0.0001,2) #0.01 可 0.001 0.0001 到0.00001 变小
proc.time()-ptm



evalbodyf5<-'
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
eval_f5<-cxxfunction(signature(Rxi='numeric',Rx='vector',Rprice='vector',Rmcs1='Matrix'),
                       body=evalbodyf5, plugin="Rcpp")


sol=
  OPTM25(0.000001,10000,In.price)$solution

sol  
eval_f5(0.01,sol,In.price,mcs1)


mcs2=mcs1-1
evalf5t<-function(x){
XP=c()
for(i in 1:1000){
  XP[i]=mcs2[,i]%*%x
}
fun=sum(exp(-xi*XP))

g=c()
for(k in 1:50){
  g[k]=sum((-xi*mcs2[k,])*exp(-xi*XP))
}
return( list( "objective"=fun, "gradient"=g ) )}




  
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
  dim(mcs1)
  
      
OPTM25=function(xi,IH,price){
    
    eval_f <- function( x ) {eval_f5(xi,x,price,mcs1)}
  
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

sol=
  OPTM25(0.01,10000,In.price)$solution

sol



eval_f5(1,3,sol,mcs1,In.price)

port.simM25=function(xi,N){
  X=matrix(0,nrow=301,ncol=N)
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=sm2[,ss]
    
    portfolio=matrix(0,nrow=5,ncol=50)
    portfolio[1,]=OPTM25(xi,10000,In.price)$solution
    X[1:61,n]=price[1:61,]%*%portfolio[1,]
    for(i in 1:4){
      IH=X[60*i+1,n]
      c.price=price[60*i+1,]
      portfolio[i+1,]=OPTM25(xi,IH,c.price)$solution
      X[(i*60+2):((i+1)*60+1),n]=price[(i*60+2):((i+1)*60+1),]%*%portfolio[i+1,]
    }
  }
  return(X)
}








# p5.z100<-port.simM25(100,500)
# p5.z10<-port.simM25(10,500)
# p5.z2<-port.simM25(2,500)
# p5.z1<-port.simM25(1,500)
# p5.5<-port.simM25(0.5,500)
# p5.1<-port.simM25(0.1,500)
p5.05<-port.simM25(0.05,500)
p5.02<-port.simM25(0.02,500)
p5.01<-port.simM25(0.01,500)
p5.005<-port.simM25(0.005,500)
p5.004<-port.simM25(0.004,500)
p5.003<-port.simM25(0.003,500)
p5.002<-port.simM25(0.002,500)
p5.001<-port.simM25(0.001,500)
p5.0005<-port.simM25(0.0005,500)
p5.0004<-port.simM25(0.0004,500)
p5.0003<-port.simM25(0.0003,500)
p5.0002<-port.simM25(0.0002,500)
p5.0001<-port.simM25(0.0001,500)
p5.00009<-port.simM25(0.00009,500)
p5.00008<-port.simM25(0.00008,500)

# f5.z100<-perfYY(p5.z100,500)
# f5.z10<-perfYY(p5.z10,500)
# f5.z2<-perfYY(p5.z2,500)
# f5.z1<-perfYY(p5.z1,500)
# f5.5<-perfYY(p5.5,500)
# f5.1<-perfYY(p5.1,500)
f5.05<-perfYY(p5.05,500)
f5.02<-perfYY(p5.02,500)
f5.01<-perfYY(p5.01,500)
f5.005<-perfYY(p5.005,500)
f5.004<-perfYY(p5.004,500)
f5.003<-perfYY(p5.003,500)
f5.002<-perfYY(p5.002,500)
f5.001<-perfYY(p5.001,500)
f5.0005<-perfYY(p5.0005,500)
f5.0004<-perfYY(p5.0004,500)
f5.0003<-perfYY(p5.0003,500)
f5.0002<-perfYY(p5.0002,500)
f5.0001<-perfYY(p5.0001,500)
f5.00009<-perfYY(p5.00009,500)
f5.00008<-perfYY(p5.00008,500)

rbind(f5.05,f5.02,f5.01,f5.005,f5.004,f5.003,f5.002,f5.001,f5.0005,f5.0004,f5.0003,f5.0002,f5.0001,f5.00009,f5.00008)
plot(density(p5.05[301,]),type='l')
lines(density(p5.001[301,]))
