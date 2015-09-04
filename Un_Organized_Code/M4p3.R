port.simM3(100,0,1000)
n=1000
MC=MC1
port.simM3=function(gamma,lambda,N){
  X=matrix(0,nrow=481,ncol=N)
  X[1,]=rep(1,N)
  portfolio1=OPTM3(10,0,1,rep(1,20))$solution
  weight=portfolio1*rep(1,20)/1
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=HAHA[,ss]
    X[2:13,n]=price[2:13,]%*%portfolio1
    for(i in 1:39){
      IH=X[12*i+1,n]
      c.price=price[12*i+1,]
      portfolio=weight*IH/c.price
      X[(i*12+2):((i+1)*12+1),n]=price[(i*12+2):((i+1)*12+1),]%*%portfolio
    }
  }
  return(X)
}

port.simM25=function(xi,N){
  X=matrix(0,nrow=481,ncol=N)
  X[1,]=rep(1,N)
  portfolio1=OPTM25(xi,1,rep(1,20))$solution
  weight=portfolio1*rep(1,20)/1
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=HAHA[,ss]
    X[2:13,n]=price[2:13,]%*%portfolio1
    for(i in 1:39){
      IH=X[12*i+1,n]
      c.price=price[12*i+1,]
      portfolio=weight*IH/c.price
      X[(i*12+2):((i+1)*12+1),n]=price[(i*12+2):((i+1)*12+1),]%*%portfolio
    }
  }
  return(X)
}
port.simM25=function(xi,N){
  X=matrix(0,nrow=481,ncol=N)
  X[1,]=rep(1,N)
  portfolio1=OPTM25(xi,1,rep(1,20))$solution
  weight=portfolio1*rep(1,20)/1
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=HAHA[,ss]
    X[2:13,n]=price[2:13,]%*%portfolio1
    for(i in 1:39){
      IH=X[12*i+1,n]
      c.price=price[12*i+1,]
      portfolio=weight*IH/c.price
      X[(i*12+2):((i+1)*12+1),n]=price[(i*12+2):((i+1)*12+1),]%*%portfolio
    }
  }
  return(X)
}

plot(density(port.simM3(25,0,n)[13,]),xlim=c(0.95,1.10),col=1)
lines(density(port.simM3(20,0,n)[13,]),col=2)
lines(density(port.simM3(15,0,n)[13,]),col=3)
lines(density(port.simM3(10,0,n)[13,]))
lines(density(port.simM3(8,0,n)[13,]),col=4)
lines(density(port.simM3(5,0,n)[13,]),col=5)
lines(density(port.simM3(4,0,n)[13,]),col=6)
lines(density(port.simM3(2,0,n)[13,]),col=1)
lines(density(port.simM3(1,0,n)[13,]))
lines(density(port.simM3(0.1,0,1000)[13,]))
lines(density(port.simM3(0.01,0,1000)[37,]))
lines(density(port.simM3(0.001,0,1000)[37,]))
lines(density(port.simM3(0.00001,0,1000)[37,]))

port.simM25=function(xi,N){
  X=matrix(0,nrow=13,ncol=N)
  X[1,]=rep(1,N)
  portfolio=OPTM25(xi,1,rep(1,20))$solution
  for(n in 1:N){
    
    ss=Css(n)#C++ 
    price=DF2[,ss]
    
    X[2:13,n]=price%*%portfolio
  }
  return(X)
}



plot(density(port.simM25(50,1000)[13,]),xlim=c(0.95,1.10),col=2)
lines(density(port.simM25(40,1000)[13,]),col=3)
lines(density(port.simM25(30,1000)[13,]),col=4)
lines(density(port.simM25(20,1000)[13,]),col=5)
lines(density(port.simM25(15,1000)[13,]),col=2)
lines(density(port.simM25(10,1000)[13,]),col=2)
lines(density(port.simM25(5,1000)[13,]),col=2)
lines(density(port.simM25(3,1000)[13,]),col=4)
MC=MC2

OI=matrix(0,nrow=481,ncol=3)
KJK=port.simM25(20,1000)
for(i in 1:481){
  OI[i,]=as.numeric(quantile(KJK[i,],probs=c(0.25,0.5,0.75)))
}
plot(0:480,OI[,1],type='l',ylim=c(0.9,2))
lines(0:480,OI[,2],col=2)
lines(0:480,OI[,3],col=3)


OI=matrix(0,nrow=481,ncol=3)
KJK=port.simM3(10,0,1000)
for(i in 1:481){
  OI[i,]=as.numeric(quantile(KJK[i,],probs=c(0.25,0.5,0.75)))
}
plot(0:480,OI[,1],type='l',ylim=c(0.9,2))
lines(0:480,OI[,2],col=2)
lines(0:480,OI[,3],col=3)
KJK


N=1000
mv50<-perfYY(port.simM3(50,0,1000),N)
mv30<-perfYY(port.simM3(30,0,1000),N)
mv25<-perfYY(port.simM3(25,0,1000),N)
mv20<-perfYY(port.simM3(20,0,1000),N)
mv17<-perfYY(port.simM3(17,0,1000),N)
mv14<-perfYY(port.simM3(14,0,1000),N)
mv12<-perfYY(port.simM3(12,0,1000),N)
mv10<-perfYY(port.simM3(10,0,1000),N)
mv8<-perfYY(port.simM3(8,0,1000),N)
mv6<-perfYY(port.simM3(6,0,1000),N)
mv4<-perfYY(port.simM3(4,0,1000),N)
mv2<-perfYY(port.simM3(2,0,1000),N)
mv1<-perfYY(port.simM3(1,0,1000),N)
mv0.9<-perfYY(port.simM3(0.9,0,1000),N)
mv0.8<-perfYY(port.simM3(0.8,0,1000),N)
M=rbind(mv30,mv25,mv20,mv17,mv14,mv12,mv10,mv8,mv6,mv4,mv2,mv1,mv0.9)
M
fs70<-perfYY(port.simM25(70,1000),N)
fs30<-perfYY(port.simM25(30,1000),N)
fs25<-perfYY(port.simM25(25,1000),N)
fs20<-perfYY(port.simM25(20,1000),N)
fs17<-perfYY(port.simM25(17,1000),N)
fs14<-perfYY(port.simM25(14,1000),N)
fs12<-perfYY(port.simM25(12,1000),N)
fs10<-perfYY(port.simM25(10,1000),N)
fs8<-perfYY(port.simM25(8,1000),N)
fs6<-perfYY(port.simM25(6,1000),N)
fs4<-perfYY(port.simM25(4,1000),N)
fs2<-perfYY(port.simM25(2,1000),N)

S=rbind(fs70,fs30,fs25,fs20,fs17,fs14,fs12,fs10,fs8,fs6)


S=MVK10[c(2,3,4,5,6,7,8,9,10),]
SSS=matrix(0,nrow=9,ncol=8)
SSS[,1]=c(30,25,20,17,14,12,10,8,6)
SSS[,2]=round(as.numeric(S[,8])*100,2)
SSS[,3]=round(as.numeric(S[,9])*100,2)
SSS[,4]=round(as.numeric(S[,1]),3)
SSS[,5]=round(as.numeric(S[,3]),3)
SSS[,6]=round(as.numeric(S[,4]),3)
SSS[,7]=round(as.numeric(S[,5]),1)
SSS[,8]=round(as.numeric(S[,7])*100,1)

SSS


mv6.10<-perfYY(port.simM3(6,10,1000),N)
mv6.50<-perfYY(port.simM3(6,50,1000),N)
mv6.100<-perfYY(port.simM3(6,100,1000),N)
mv6.500<-perfYY(port.simM3(6,500,1000),N)
mv6.1000<-perfYY(port.simM3(6,1000,1000),N)
mv6.2000<-perfYY(port.simM3(6,2000,1000),N)
mv6.3000<-perfYY(port.simM3(6,3000,1000),N)
mv6.4000<-perfYY(port.simM3(6,4000,1000),N)
mv6.5000<-perfYY(port.simM3(6,5000,1000),N)
mv6.10000<-perfYY(port.simM3(6,10000,1000),N)
mv6.20000<-perfYY(port.simM3(6,20000,1000),N)
mvk6<-rbind(mv6,mv6.10,mv6.100,mv6.1000,mv6.2000,mv6.3000,mv6.4000,mv6.5000,mv6.10000,mv6.20000)

mvk6


m50<-perfYY(port.simM3(50,10,1000),N)
m30<-perfYY(port.simM3(30,10,1000),N)
m25<-perfYY(port.simM3(25,10,1000),N)
m20<-perfYY(port.simM3(20,10,1000),N)
m17<-perfYY(port.simM3(17,10,1000),N)
m14<-perfYY(port.simM3(14,10,1000),N)
m12<-perfYY(port.simM3(12,10,1000),N)
m10<-perfYY(port.simM3(10,10,1000),N)
m8<-perfYY(port.simM3(8,10,1000),N)
m6<-perfYY(port.simM3(6,10,1000),N)
m4<-perfYY(port.simM3(4,10,1000),N)
m2<-perfYY(port.simM3(2,10,1000),N)
m1<-perfYY(port.simM3(1,10,1000),N)

MVK10=rbind(m30,m25,m20,m17,m14,m12,m10,m8,m6,m4,m2,m1)



m150<-perfYY(port.simM3(50,1,1000),N)
m130<-perfYY(port.simM3(30,1,1000),N)
m125<-perfYY(port.simM3(25,1,1000),N)
m120<-perfYY(port.simM3(20,1,1000),N)
m117<-perfYY(port.simM3(17,1,1000),N)
m114<-perfYY(port.simM3(14,1,1000),N)
m112<-perfYY(port.simM3(12,1,1000),N)
m110<-perfYY(port.simM3(10,1,1000),N)
m18<-perfYY(port.simM3(8,1,1000),N)
m16<-perfYY(port.simM3(6,1,1000),N)
m14<-perfYY(port.simM3(4,1,1000),N)
m12<-perfYY(port.simM3(2,1,1000),N)
m11<-perfYY(port.simM3(1,1,1000),N)

MVK1=rbind(m150,m130,m125,m120,m117,m114,m112,m110,m18,m16,m14,m12,m11)
MVK1

m150<-perfYY(port.simM3(50,100,1000),N)
m130<-perfYY(port.simM3(30,100,1000),N)
m125<-perfYY(port.simM3(25,100,1000),N)
m120<-perfYY(port.simM3(20,100,1000),N)
m117<-perfYY(port.simM3(17,100,1000),N)
m114<-perfYY(port.simM3(14,100,1000),N)
m112<-perfYY(port.simM3(12,100,1000),N)
m110<-perfYY(port.simM3(10,100,1000),N)
m18<-perfYY(port.simM3(8,100,1000),N)
m16<-perfYY(port.simM3(6,100,1000),N)
m14<-perfYY(port.simM3(4,100,1000),N)
m12<-perfYY(port.simM3(2,100,1000),N)
m11<-perfYY(port.simM3(1,100,1000),N)

MVK1=rbind(m130,m125,m120,m117,m114,m112,m110,m18,m16,m14,m12,m11)
MVK1




plot(S[,3],S[,1],type='l',ylim=c(0,0.18))
lines(M[,3],M[,1],col=2)
lines(MVK1[,3],MVK1[,1],col=3)
lines(MVK10[,3],MVK10[,1],col=4)



oooo=function(seed,N){
  MCE=c();MCVA=c();MCK=c()
  for(i in 1:20){
    sim=ugarchsim(res[[i]], n.sim = 24, m.sim = N, n.start=100,rseed = seed)
    JK=colSums(fitted(sim)[13:24,])
    M=mean(exp(JK))
    while(M>1.2){
      sim=ugarchsim(res[[i]], n.sim = 24, m.sim = N, n.start=100,rseed = seed)
      JK=colSums(fitted(sim)[13:24,])
      M=mean(exp(JK))
    }
    MCE[i]=mean(exp(JK))
    MCVA[i]=var(exp(JK))
    MCK[i]=mean((exp(JK)-MCE[i])^4)
  }
  return(list('E'=MCE,'VA'=MCVA,'K'=MCK))
}


ZX=oooo(NA,10000)
S
plot(0:12,port.simM25(50,1000)[,12],type='l',ylim=c(0.96,1.02))
abline(h=1,lty=2)
lines(0:12,port.simM25(40,1000)[,12],col=2)
lines(0:12,port.simM25(30,1000)[,12],col=3)
lines(0:12,port.simM25(20,1000)[,12],col=4)
lines(0:12,port.simM25(10,1000)[,12],col=5)
lines(0:12,port.simM25(6,1000)[,12],col=6)

HAHA=MKM()
plot(0:200,c(1,HAHA%*%OPTM3(30,0,1,rep(1,20))$solution)-1,type='l',ylim=c(-0.15,0.15),main='Mean-Vairance',ylab='profit',xlab='time')
abline(h=0,lty=2)
lines(0:200,c(1,HAHA%*%OPTM3(25,0,1,rep(1,20))$solution)-1,col=2)
lines(0:200,c(1,HAHA%*%OPTM3(20,0,1,rep(1,20))$solution)-1,col=3)
lines(0:200,c(1,HAHA%*%OPTM3(15,0,1,rep(1,20))$solution)-1,col=4)
lines(0:200,c(1,HAHA%*%OPTM3(10,0,1,rep(1,20))$solution)-1,col=5)
lines(0:200,c(1,HAHA%*%OPTM3(5,0,1,rep(1,20))$solution)-1,col=6)



plot(0:200,c(1,HAHA%*%OPTM25(60,1,rep(1,20))$solution)-1,type='l',ylim=c(-0.15,0.15),main='Risk Sensitive',ylab='profit',xlab='time')
abline(h=0,lty=2)
lines(0:200,c(1,HAHA%*%OPTM25(50,1,rep(1,20))$solution)-1,col=2)
lines(0:200,c(1,HAHA%*%OPTM25(40,1,rep(1,20))$solution)-1,col=3)
lines(0:200,c(1,HAHA%*%OPTM25(30,1,rep(1,20))$solution)-1,col=4)
lines(0:200,c(1,HAHA%*%OPTM25(20,1,rep(1,20))$solution)-1,col=5)
lines(0:200,c(1,HAHA%*%OPTM25(10,1,rep(1,20))$solution)-1,col=6)

length(MVK1[,3])
x1<-as.numeric(c(M[,9],S[,9],MVK10[,9],MVK1[,9]))
x2<-as.numeric(c(M[,8],S[,8],MVK10[,8],MVK1[,8]))
x3=c(rep('Mean-Variance',13),rep('Risk-Sensitive',10),rep('MVK fixLambda 10',12),rep('MVK fixLambda 100',12))

k<-data.frame(x1,x2,x3)
colnames(k) <- c("Maximum_Drawdown", "Sharpe_Ratio","Cost_function")
k
p1<-ggplot(data=k,aes(x=Maximum_Drawdown,y=Sharpe_Ratio,color=factor(Cost_function)))+geom_point()+ stat_smooth()
p1
k
plot(density(rt(10000,5)))
plot(density(rexp(100000,1)))
