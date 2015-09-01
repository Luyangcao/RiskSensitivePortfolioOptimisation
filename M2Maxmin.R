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

################################

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



#####

# lb <- -(0.25*1)/abs(In.price)
# ub <- (0.25*1)/abs(In.price)
# x0 <- runif(19,lb[1:19],ub[1:19])
# 
# cons=function(x0){
# x20=(1-x0%*%In.price[1:19])/In.price[20]
# xx=c(x0,x20)
# return(xx)}
# 
# rnorm(19,0,1)
# func=function(xx){
# MaMi(20,0,0,xx,MCm2$E,MCm2$VA,MCm2$K,mxx,In.price)
# }
# 
# 
# XX=cons(runif(19,lb[1:19],ub[1:19]))
# 
# 
# simulated_annealing <- function(func, s0, niter, step ) {
#   
#   # Initialize
#   ## s stands for state
#   ## f stands for function value
#   ## b stands for best
#   ## c stands for current
#   ## n stands for neighbor
# 
#   s_b <- s_c <- s_n <- s0
#   f_b <- f_c <- f_n <- func(s_n)
#   message("It\tBest\tCurrent\tNeigh\tTemp")
#   message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", 0L, f_b, f_c, f_n, 1))
#   
#   for (k in 1:niter) {     
#     Temp <- (1 - step)^k
#     # consider a random neighbor
#     s_n <- runif(x,lb,ub)
#     f_n <- func(s_n)
#     # update current state
#     if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
#       s_c <- s_n
#       f_c <- f_n
#     }
#     # update best state
#     if (f_n < f_b) {
#       s_b <- s_n
#       f_b <- f_n         
#     }
#     message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", k, f_b, f_c, f_n, Temp))
#   }
#   return(list(iterations = niter, best_value = f_b, best_state = s_b))
# }
# XX=cons(runif(19,lb[1:19],ub[1:19]))
# XX=bound(XX)
# LL=simulated_annealing(func,XX,niter=1000,step=0.00001)
# LL
# 
# bound=function(XX){
# OX=cons(rnorm(19,XX,0.001))
# while(sum(OX>ub)+sum(OX<lb)>0){
#   OX=cons(rnorm(19,XX,0.001))
# }
# return(OX)}
# 
# XX
# OX=bound(XX)
# 
# OPTM2(20,0,1,In.price)$solution
# LL$best_state
# 
# 

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



LL=OP2(25,0,2,1,In.price)
LL$solution%*%In.price
LL
port.simM23=function(gamma,lambda,theta,N){
  X=matrix(0,nrow=91,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=20)
  portfolio[1,]=OP2(gamma,lambda,theta,1,In.price)$solution
  for(n in 1:N){
    ss=Css(n)#C++ 
    price=sm2W[,ss]
    
    
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
