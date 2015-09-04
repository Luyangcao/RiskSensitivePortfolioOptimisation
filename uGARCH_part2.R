
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






## mcs or mcs1 is simulated terminal value of each stock, used to monte carlo method in optimization

SSIM=function(){
  SIM=matrix(0,nrow=12,ncol=20000)
  for(i in 1:20){
    sim=ugarchsim(res[[i]], n.sim = 12, m.sim = 1000, n.start=100)
    JK=matrix(0,nrow=12,ncol=1000)
    for(j in 1:1000){
      JK[,j]=cumsum(fitted(sim)[,j])}
    SIM[,((i-1)*1000+1):(i*1000)]=exp(JK)
  }
  return(SIM)
}
DF=SSIM()
kk=DF[12,]
kkk=matrix(0,nrow=20,ncol=1000)
for(i in 1:1000){
  cl=Css(i)
  kkk[,i]=kk[cl]
}
tmcs1=kkk




##########################################################################################################

### risk sensitive portfolio optmizer. 
# xi is risk tolerance parameter , IH=initial holding
OPTM25=function(xi,IH,price){
  
  eval_f <- function( x ) {teval_f5(xi,x,price,tmcs1)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.3*IH)/abs(price)
  ub <- (0.3*IH)/abs(price)
  set.seed(998)
  x0 <- runif(20,lb,ub)
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-8 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-8,
                "maxeval" = 10000,
                "local_opts" = local_opts )
  res0 <- nloptr(x0,eval_f=eval_f,lb=lb,ub=ub,eval_g_ineq=eval_g_eq,opts=opts)
  return(res0)}







###############################################################################################################

#### calculate portfolio value of GBM model of risk sensitive reward function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  port.simM25=function(xi,N){
    X=matrix(0,nrow=13,ncol=N)
    X[1,]=rep(1,N)
    portfolio=OPTM25(xi,1,rep(1,40))$solution
    for(n in 1:N){
      
      ss=Css(n)#C++ 
      price=DF2[,ss]
      
      X[2:13,n]=price%*%portfolio
    }
    return(X)
  }


