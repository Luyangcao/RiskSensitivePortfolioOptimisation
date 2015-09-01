
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

Maxsim<-tmcsim(best,rep(1,50))




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
NumericMatrix Vs(37,50);


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
NumericVector XP(37);
for(int i=1; i<=37; i++){
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





OP2=function(gamma,lambda,theta,IH,price){
  E=tMCm2$E*price
  VA=tMCm2$VA*(price^2)
  K=tMCm2$K*(price^4)
  
  
  eval_f <- function(xx){
    MaMi(gamma,lambda,theta,xx,E,VA,K,Maxsim,price)
  }
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    return( constr)
  }
  
  # initial values and bounds
  lb <- -(0.2*IH)/abs(price)
  ub <- (0.2*IH)/abs(price)
  
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



LL=OP2(25,0,2,1,rep(1,50))
LL$solution%*%rep(1,50)
LL

port.Max=function(gamma,lambda,theta,N){
  X=matrix(0,nrow=37,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=50)
  portfolio[1,]=OP2(gamma,lambda,theta,1,rep(1,50))$solution
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