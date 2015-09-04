
best=list()

Black <- read.csv("~/Desktop/Project/prices/msc_thesis_data.csv")
summary(Black$asset_id)


riskfree<-c(1.51,1.57,1.65,1.72,1.62,1.78,2.02,1.89,1.87,1.70,2.16,2.52,2.60,2.78,2.64,2.57,2.75,3.04,
            2.67,2.66,2.73,2.67,2.48,2.53,2.58,2.35,2.52,2.35,2.18,2.17,
            1.68,2.00,1.94,2.05,2.12,2.35)
RF=rep(0,50)
z <- as.yearmon(1960 + seq(0, 11)/12)
RF=riskfree/1200
Rf=data.frame(RF,row.names=as.Date(z))
mean(RF)

perfYY=function(YYY,N){
  X=YYY
  SR=numeric(N);IR=numeric(N);ADD=numeric(N);MDD=numeric(N);NDD=numeric(N);MT=numeric(N);LS=numeric(N);Re=numeric(N)
  for(i in 1:N){
    k=X[,i]
    ret=diff(k)/k[1:36] #return
    R=data.frame(ret,row.names=as.Date(z))
    SR[i]=SharpeRatio(R,Rf,FUN="StdDev") # 
    ADD[i]=AverageDrawdown(R) #
    MDD[i]=maxDrawdown(R)#
    FD=findDrawdowns(R)
    NDD[i]=length(FD$from)
    MT[i]=max(FD$to-FD$from)
    LS[i]=as.numeric(k[37]<1)
    Re[i]=k[37]-1
    IR[i]=InformationRatio(R,Rf)
  }
  return(list('SR'=mean(SR),'IR'=mean(IR),'MDD'=mean(MDD),'ADD'=mean(ADD),'NDD'=mean(NDD),
              "MT"=mean(MT),"LS"=mean(LS),"R"=mean(Re),"V"=sd(Re)))}






A102=Black[which(Black$asset_id=='Asset_102'),]
A105=Black[which(Black$asset_id=='Asset_105'),]
A109=Black[which(Black$asset_id=='Asset_109'),]
A116=Black[which(Black$asset_id=='Asset_116'),]
A118=Black[which(Black$asset_id=='Asset_118'),]
A122=Black[which(Black$asset_id=='Asset_122'),]
A127=Black[which(Black$asset_id=='Asset_127'),]
A132=Black[which(Black$asset_id=='Asset_132'),]
A134=Black[which(Black$asset_id=='Asset_134'),]
A136=Black[which(Black$asset_id=='Asset_136'),]
A138=Black[which(Black$asset_id=='Asset_138'),]
A14=Black[which(Black$asset_id=='Asset_14'),]
A141=Black[which(Black$asset_id=='Asset_141'),]
A142=Black[which(Black$asset_id=='Asset_142'),]
A150=Black[which(Black$asset_id=='Asset_150'),]
A158=Black[which(Black$asset_id=='Asset_158'),]
A160=Black[which(Black$asset_id=='Asset_160'),]
A161=Black[which(Black$asset_id=='Asset_161'),]
A165=Black[which(Black$asset_id=='Asset_165'),]
A175=Black[which(Black$asset_id=='Asset_175'),]
A177=Black[which(Black$asset_id=='Asset_177'),]
A60=Black[which(Black$asset_id=='Asset_60'),]
A185=Black[which(Black$asset_id=='Asset_185'),]
A19=Black[which(Black$asset_id=='Asset_19'),]
A20=Black[which(Black$asset_id=='Asset_20'),]
A22=Black[which(Black$asset_id=='Asset_22'),]
A29=Black[which(Black$asset_id=='Asset_29'),]
A30=Black[which(Black$asset_id=='Asset_30'),]
A32=Black[which(Black$asset_id=='Asset_32'),]
A33=Black[which(Black$asset_id=='Asset_33'),]
A35=Black[which(Black$asset_id=='Asset_35'),]
A4=Black[which(Black$asset_id=='Asset_4'),]
A40=Black[which(Black$asset_id=='Asset_40'),]
A43=Black[which(Black$asset_id=='Asset_43'),]
A5=Black[which(Black$asset_id=='Asset_5'),]
A52=Black[which(Black$asset_id=='Asset_52'),]
A54=Black[which(Black$asset_id=='Asset_54'),]
A57=Black[which(Black$asset_id=='Asset_57'),]
A58=Black[which(Black$asset_id=='Asset_58'),]
A61=Black[which(Black$asset_id=='Asset_61'),]
A64=Black[which(Black$asset_id=='Asset_64'),]
A72=Black[which(Black$asset_id=='Asset_72'),]
A8=Black[which(Black$asset_id=='Asset_8'),]
A80=Black[which(Black$asset_id=='Asset_80'),]
A83=Black[which(Black$asset_id=='Asset_83'),]
A84=Black[which(Black$asset_id=='Asset_84'),]
A85=Black[which(Black$asset_id=='Asset_85'),]
A93=Black[which(Black$asset_id=='Asset_93'),]
A99=Black[which(Black$asset_id=='Asset_99'),]
A77=Black[which(Black$asset_id=='Asset_77'),]


M4=cbind(A102$realised_signal,A105$realised_signal,A109$realised_signal,A116$realised_signal
             ,A118$realised_signal,A122$realised_signal,A127$realised_signal,
      A132$realised_signal,A134$realised_signal,A136$realised_signal,A138$realised_signal,
      A14$realised_signal,A141$realised_signal,A142$realised_signal,A150$realised_signal,
      A158$realised_signal,A160$realised_signal,A161$realised_signal,A165$realised_signal,
      A175$realised_signal,A177$realised_signal,A60$realised_signal,A185$realised_signal,A19$realised_signal,A20$realised_signal,
      A22$realised_signal,A29$realised_signal,A30$realised_signal,A32$realised_signal,A33$realised_signal,
      A35$realised_signal,A4$realised_signal,A40$realised_signal,A43$realised_signal,A5$realised_signal,
      A52$realised_signal,A54$realised_signal,A57$realised_signal,A58$realised_signal,A61$realised_signal,
      A64$realised_signal,A72$realised_signal,A8$realised_signal,A80$realised_signal,A83$realised_signal,
      A84$realised_signal,A85$realised_signal,A93$realised_signal,A99$realised_signal,A77$realised_signal)


BBlack=cbind(A102$spec_signal,A105$spec_signal,A109$spec_signal,A116$spec_signal
             ,A118$spec_signal,A122$spec_signal,A127$spec_signal,
             A132$spec_signal,A134$spec_signal,A136$spec_signal,A138$spec_signal,
             A14$spec_signal,A141$spec_signal,A142$spec_signal,A150$spec_signal,
             A158$spec_signal,A160$spec_signal,A161$spec_signal,A165$spec_signal,
             A175$spec_signal,A177$spec_signal,A60$spec_signal,A185$spec_signal,A19$spec_signal,A20$spec_signal,
             A22$spec_signal,A29$spec_signal,A30$spec_signal,A32$spec_signal,A33$spec_signal,
             A35$spec_signal,A4$spec_signal,A40$spec_signal,A43$spec_signal,A5$spec_signal,
             A52$spec_signal,A54$spec_signal,A57$spec_signal,A58$spec_signal,A61$spec_signal,
             A64$spec_signal,A72$spec_signal,A8$spec_signal,A80$spec_signal,A83$spec_signal,
             A84$spec_signal,A85$spec_signal,A93$spec_signal,A99$spec_signal,A77$spec_signal)

BBlack

  for(i in 1:50){
    y=BBlack[1:72,1]
    rett<-log(y+1)
    fit <- stepAIC.ghyp(rett, dist = c( "hyp", "NIG", "t", "gauss"),symmetric = NULL, control = list(maxit = 5000),na.rm = T,silent = TRUE)
    bestfit <-fit$best.model
    best[i]=bestfit}
best






tcxxsm2<-'
Environment ghyp("package:ghyp");
Function rghyp = ghyp["rghyp"];
NumericMatrix sm2(37,500000);

List best(Rbest);
NumericVector  price(Rprice);

for (int i=1; i<=50; i++) {
for(int m=1; m<=10000; m++){
sm2(1-1,(i-1)*10000+m-1)=price[i-1];
NumericVector  gg(36);
gg=rghyp(36,best[i-1]);
for(int n=1; n<=36; n++){
sm2(n+1-1,(i-1)*10000+m-1)=sm2(n-1,(i-1)*10000+m-1)*exp(gg[n-1]);
}
}
}
return wrap(sm2);'

tCppsm2<-cxxfunction(signature(Rbest='List',Rprice='vector'),
                    body=tcxxsm2, plugin="Rcpp")


ptm<-proc.time()
tsm2W=tCppsm2(best,rep(1,50))
proc.time()-ptm


tMcsimbodym2<-'
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
NumericMatrix pricesim(13,100000*50);



for(int i=1;i<=50;i++){
for(int m=1; m<=100000;m++){
pricesim(1-1,(i-1)*100000+m-1)=IV;
NumericVector  gg(12);
gg=rghyp(12,best[i-1]);
for(int n=1;n<=12;n++){
pricesim(n+1-1,(i-1)*100000+m-1)=pricesim(n-1,(i-1)*100000+m-1)*exp(gg[n-1]);
}
}}

NumericVector E(50);
NumericVector VA(50);
NumericVector K(50);


for(int i=1; i<=50;i++){
NumericVector st(100000);
for(int m=1; m<=100000;m++){
st[m-1]=pricesim(13-1,(i-1)*100000+m-1);
}
E[i-1]=as<double>(mean(st));
VA[i-1]=as<double>(var(st));
NumericVector st2(100000);
for(int n=1; n<=100000;n++){
st2[n-1]=(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1])*(st[n-1]-E[i-1]);}
K[i-1]=as<double>(mean(st2));
}
List ret; ret["E"] = E; ret["VA"] = VA;ret["K"] = K;
return wrap(ret);
'
twwm2<-cxxfunction(signature(RIV='numeric',Rbest='List'),
                  body=tMcsimbodym2, plugin="Rcpp")
tMCm2<-twwm2(1,best)
tMCm21<-twwm2(1,best)
tMCm22<-twwm2(1,best)
tMCm2.100<-twwm2(100,best)

tMCm21$K
tMCm2$K
tMCm2.100$K


evalbody4<-'
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
eval_f4<-cxxfunction(signature(Rgamma='numeric',Rlambda='numeric',Rx='vector',RE='vector',RVA='vector',RK='vector'),
                     body=evalbody4, plugin="Rcpp")








OPTM3(100,0,1,rep(1,30))$solution%*%rep(1,30)

OPTM3=function(gamma,lambda,IH,price){
  
  E=tMCm2$E*price
  VA=tMCm2$VA*price^2
  K=tMCm2$K*price^4
  
  eval_f <- function( x ) {eval_f4(gamma,lambda,x,E,VA,K)}
  
  eval_g_eq <- function( x ) {
    constr <- c( x%*%price -IH )
    grad <- c( price)
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  # initial values and bounds
  lb <- -(0.25*IH)/abs(price)
  ub <- (0.25*IH)/abs(price)
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

ssbody<-'
    double n = as<double>(Rn);
NumericVector ss(50);
for(int i=1; i<=50; i++){
ss(i-1)=(i-1)*10000+n;
}
return wrap(ss);
'
Css<-cxxfunction(signature(Rn='numeric'),
                 body=ssbody, plugin="Rcpp")

var(rt(1000000,4))
var(rnorm(1000,0,sqrt(2)))

port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=37,ncol=N)
  portfolio=OPTM3(gamma,lambda,1,rep(1,50))$solution
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=tsm2W[,ss]
    
    X[,n]=price%*%portfolio
  }
  return(X)
}


port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=37,ncol=N)
  portfolio=matrix(0,nrow=3,ncol=50)
  portfolio[1,]=OPTM3(gamma,lambda,1,rep(1,50))$solution
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