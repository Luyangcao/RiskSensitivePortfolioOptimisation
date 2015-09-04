#data
library(quadprog)
Stocks <- read.csv("~/Desktop/work2/Stocks")
Stocks<-Stocks[,-1]
Four<-apply(Stocks, 2, rev)

Data1<-cbind(Four[1:267,],Four[2:268,])

colnames(Data1)<-c("xs1","xs2","xs3","xs4","ys1","ys2","ys3","ys4")
Data1<-data.frame(Data1)
##########################################################################################
###estimate the parameter of geometric bronian motion
##########################################################################################
my=Data1[1:177,5:8]
mx=Data1[1:177,1:4]


mgbm=my/mx
loggbm=(log(mgbm))
sd1=sqrt(var(lm(loggbm[,1]~1)$residuals))
mu1=lm(loggbm[,1]~1)$coefficient+sd1^2/2

sd2=sqrt(var(lm(loggbm[,2]~1)$residuals))
mu2=lm(loggbm[,2]~1)$coefficient+sd2^2/2

sd3=sqrt(var(lm(loggbm[,3]~1)$residuals))
mu3=lm(loggbm[,3]~1)$coefficient+sd3^2/2

sd4=sqrt(var(lm(loggbm[,4]~1)$residuals))
mu4=lm(loggbm[,4]~1)$coefficient+sd4^2/2


####################################################
####################peroid 1
#############################################
IP=Four[178,]
##############

opt3<-function(IP,TV,C,IH,gamma){
  Amat <- cbind(matrix(c(100,IP)));
  bvec<-c(sum(IH));
  k<-solve.QP(C,0.5/gamma*TV,Amat,bvec,meq=1)
  return(k)
}


library(quadprog)



#####
##########################################################################################
###############Simulation
##########################################################################################


##########



s1=matrix(0,nrow=100000,ncol=5)
for(i in 1:100000){
  s1[i,2]=IP[1]*exp((mu1-sd1^2/2)*90+sd1*rnorm(1,0,sqrt(90)))
  s1[i,3]=IP[2]*exp((mu2-sd2^2/2)*90+sd2*rnorm(1,0,sqrt(90)))
  s1[i,4]=IP[3]*exp((mu3-sd3^2/2)*90+sd3*rnorm(1,0,sqrt(90)))
  s1[i,5]=IP[4]*exp((mu4-sd4^2/2)*90+sd4*rnorm(1,0,sqrt(90)))
  s1[i,1]=100*exp(90*0.0001)
}
Es=matrix(c(mean(s1[,1]),mean(s1[,2]),mean(s1[,3]),mean(s1[,4]),mean(s1[,5])),ncol=1)
vs=diag(c(0.000000001,var(s1[,2]),var(s1[,3]),var(s1[,4]),var(s1[,5])))
k1=s1%*%opt3(IP,Es,vs,10000,0.1)$solution
s1=matrix(0,nrow=100000,ncol=5)

for(i in 1:100000){
  s1[i,2]=IP[1]*exp((mu1-sd1^2/2)*90+sd1*rnorm(1,0,sqrt(90)))
  s1[i,3]=IP[2]*exp((mu2-sd2^2/2)*90+sd2*rnorm(1,0,sqrt(90)))
  s1[i,4]=IP[3]*exp((mu3-sd3^2/2)*90+sd3*rnorm(1,0,sqrt(90)))
  s1[i,5]=IP[4]*exp((mu4-sd4^2/2)*90+sd4*rnorm(1,0,sqrt(90)))
  s1[i,1]=100*exp(90*0.0001)
}
Es=matrix(c(mean(s1[,1]),mean(s1[,2]),mean(s1[,3]),mean(s1[,4]),mean(s1[,5])),ncol=1)
vs=diag(c(0.000000001,var(s1[,2]),var(s1[,3]),var(s1[,4]),var(s1[,5])))
k2=s1%*%opt3(IP,Es,vs,10000,0.01)$solution
s1=matrix(0,nrow=100000,ncol=5)

for(i in 1:100000){
  s1[i,2]=IP[1]*exp((mu1-sd1^2/2)*90+sd1*rnorm(1,0,sqrt(90)))
  s1[i,3]=IP[2]*exp((mu2-sd2^2/2)*90+sd2*rnorm(1,0,sqrt(90)))
  s1[i,4]=IP[3]*exp((mu3-sd3^2/2)*90+sd3*rnorm(1,0,sqrt(90)))
  s1[i,5]=IP[4]*exp((mu4-sd4^2/2)*90+sd4*rnorm(1,0,sqrt(90)))
  s1[i,1]=100*exp(90*0.0001)
}
Es=matrix(c(mean(s1[,1]),mean(s1[,2]),mean(s1[,3]),mean(s1[,4]),mean(s1[,5])),ncol=1)
vs=diag(c(0.000000001,var(s1[,2]),var(s1[,3]),var(s1[,4]),var(s1[,5])))
k3=s1%*%opt3(IP,Es,vs,10000,0.05)$solution
s1=matrix(0,nrow=100000,ncol=5)


for(i in 1:100000){
  s1[i,2]=IP[1]*exp((mu1-sd1^2/2)*90+sd1*rnorm(1,0,sqrt(90)))
  s1[i,3]=IP[2]*exp((mu2-sd2^2/2)*90+sd2*rnorm(1,0,sqrt(90)))
  s1[i,4]=IP[3]*exp((mu3-sd3^2/2)*90+sd3*rnorm(1,0,sqrt(90)))
  s1[i,5]=IP[4]*exp((mu4-sd4^2/2)*90+sd4*rnorm(1,0,sqrt(90)))
  s1[i,1]=100*exp(90*0.0001)
}
Es=matrix(c(mean(s1[,1]),mean(s1[,2]),mean(s1[,3]),mean(s1[,4]),mean(s1[,5])),ncol=1)
vs=diag(c(0.000000001,var(s1[,2]),var(s1[,3]),var(s1[,4]),var(s1[,5])))
k4=s1%*%opt3(IP,Es,vs,10000,0.005)$solution
###########################################################################
####plot
###########################################################################
###########################################################################
plot(density(k1),xlab='',ylab='',main="",col=1,xlim=c(9900,10800))
abline(v=mean(k1),col=1,lty=2)
lines(density(k2),col=2)
abline(v=mean(k2),col=2,lty=2)
lines(density(k3),col="green")
abline(v=mean(k3),col=3,lty=2)
lines(density(k4),col=4)
abline(v=mean(k4),col=4,lty=2)
legend('topright', c("0.1","0.5","0.05","0.005") , 
       lty=1, col=c("black",'red', 'green', 'blue'), bty='n', cex=1)


####################################################
####################peroid~
#############################################
p1<-opt3(IP,Es,vs,10000,0.1)$solution
p2<-opt3(IP,Es,vs,10000,0.01)$solution
p3<-opt3(IP,Es,vs,10000,0.05)$solution
p4<-opt3(IP,Es,vs,10000,0.005)$solution
###########
#############

IP=Four[268,]
##############
zIH1<-p1%*%c(100.9041,IP)
zIH2<-p2%*%c(100.9041,IP)
zIH3<-p3%*%c(100.9041,IP)
zIH4<-p4%*%c(100.9041,IP)

##
zIH=c(zIH1,zIH2,zIH3,zIH4)

zIH
