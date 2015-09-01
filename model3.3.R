N=1000
mv100<-perfYY(port.simM3(100,0,1000),N)
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

M=rbind(mv100,mv30,mv25,mv20,mv17,mv14,mv12,mv10,mv8,mv6,mv4,mv2,mv1)
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


S=rbind(fs70,fs30,fs25,fs20,fs17,fs14,fs12,fs10,fs8,fs6,fs4,fs2)
S
plot(S[,3],S[,1],type='l',ylim=c(0,0.5),xlim=c(0,0.08))
lines(M[,3],M[,1],col=2)

lines(K[,3],K[,1],col=4)
M[,c(1,3)]
S[,c(1,3)]




mv6.10<-perfYY(port.simM3(6,10,1000),1000)
mv6.50<-perfYY(port.simM3(6,50,1000),1000)
mv6.200<-perfYY(port.simM3(6,200,1000),1000)
mv6.500<-perfYY(port.simM3(6,500,1000),1000)
mv6.1000<-perfYY(port.simM3(6,1000,1000),1000)
mv6.2000<-perfYY(port.simM3(6,2000,1000),1000)
mv6.3000<-perfYY(port.simM3(6,3000,1000),1000)
mv6.4000<-perfYY(port.simM3(6,4000,1000),1000)
mv6.5000<-perfYY(port.simM3(6,5000,1000),1000)
mv6.10000<-perfYY(port.simM3(6,10000,1000),1000)
mv6.20000<-perfYY(port.simM3(6,20000,1000),1000)


K=rbind(mv6.10,mv6.50,mv6.200,mv6.500,mv6.1000,mv6.2000,mv6.3000,mv6.4000,mv6.5000,mv6.10000,mv6.20000)
K




mv4.10<-perfYY(port.simM3(4,10,1000),1000)
mv4.50<-perfYY(port.simM3(4,50,1000),1000)
mv4.200<-perfYY(port.simM3(4,200,1000),1000)
mv4.500<-perfYY(port.simM3(4,500,1000),1000)
mv4.1000<-perfYY(port.simM3(4,1000,1000),1000)
mv4.2000<-perfYY(port.simM3(4,2000,1000),1000)
mv4.3000<-perfYY(port.simM3(4,3000,1000),1000)
mv4.4000<-perfYY(port.simM3(4,4000,1000),1000)
mv4.5000<-perfYY(port.simM3(4,5000,1000),1000)
mv4.10000<-perfYY(port.simM3(4,10000,1000),1000)
mv4.20000<-perfYY(port.simM3(4,20000,1000),1000)



K4=rbind(mv4.10,mv4.50,mv4.200,mv4.500,mv4.1000,mv4.2000,mv4.3000,mv4.4000,mv4.5000,mv4.10000,mv4.20000)




mv10.10<-perfYY(port.simM3(10,10,1000),1000)
mv10.50<-perfYY(port.simM3(10,50,1000),1000)
mv10.200<-perfYY(port.simM3(10,200,1000),1000)
mv10.500<-perfYY(port.simM3(10,500,1000),1000)
mv10.1000<-perfYY(port.simM3(10,1000,1000),1000)
mv10.2000<-perfYY(port.simM3(10,2000,1000),1000)
mv10.3000<-perfYY(port.simM3(10,3000,1000),1000)
mv10.4000<-perfYY(port.simM3(10,4000,1000),1000)
mv10.5000<-perfYY(port.simM3(10,5000,1000),1000)
mv10.10000<-perfYY(port.simM3(10,10000,1000),1000)
mv10.20000<-perfYY(port.simM3(10,20000,1000),1000)



K10=rbind(mv10,mv10.10,mv10.200,mv10.500,mv10.1000,mv10.2000,mv10.3000,mv10.5000,mv10.10000,mv10.20000)
K10



max4.f1<-perfYY(port.Max(4,0,-1,1000),1000)
max4.0.01<-perfYY(port.Max(4,0,0.01,1000),1000)
max4.0.1<-perfYY(port.Max(4,0,0.1,1000),1000)
max4.1<-perfYY(port.Max(4,0,1,1000),1000)
max4.2<-perfYY(port.Max(4,0,2,1000),1000)
max4.4<-perfYY(port.Max(4,0,4,1000),1000)
max4.6<-perfYY(port.Max(4,0,6,1000),1000)
max4.8<-perfYY(port.Max(4,0,8,1000),1000)
max4.10<-perfYY(port.Max(4,0,10,1000),1000)
max4.12<-perfYY(port.Max(4,0,12,1000),1000)
max4.15<-perfYY(port.Max(4,0,15,1000),1000)


Max4=rbind(mv4,max4.0.1,max4.1,max4.2,max4.4,max4.6,max4.8,max4.10,max4.12,max4.15)
Max4

N=1000
mxv100<-perfYY(port.Max(100,0,0.1,1000),N)
mxv30<-perfYY(port.Max(30,0,0.1,1000),N)
mxv25<-perfYY(port.Max(25,0,0.1,1000),N)
mxv20<-perfYY(port.Max(20,0,0.1,1000),N)
mxv17<-perfYY(port.Max(17,0,0.1,1000),N)
mxv14<-perfYY(port.Max(14,0,0.1,1000),N)
mxv12<-perfYY(port.Max(12,0,0.1,1000),N)
mxv10<-perfYY(port.Max(10,0,0.1,1000),N)
mxv8<-perfYY(port.Max(8,0,0.1,1000),N)
mxv6<-perfYY(port.Max(6,0,0.1,1000),N)
mxv4<-perfYY(port.Max(4,0,0.1,1000),N)
mxv2<-perfYY(port.Max(2,0,0.1,1000),N)
mxv1<-perfYY(port.Max(1,0,0.1,1000),N)

MXV=rbind(mxv100,mxv30,mxv25,mxv20,mxv17,mxv14,mxv12,mxv10,mxv8,mxv6,mxv4,mxv2,mxv1)
MXV
length(MXV[,1])
x1<-as.numeric(c(M[,3],S[,3],K[,3],K4[,3],K10[,3],Max4[,3],MXV[,3]))
x2<-as.numeric(c(M[,1],S[,1],K[,1],K4[,1],K10[,1],Max4[,1],MXV[,1]))
x3=c(rep('Mean-Variance',13),rep('Risk-Sensitive',12),rep('MVK6',11),rep('MVK4',11),rep('MVK10',10),rep('Max4',10),rep('MVX',13))

k<-data.frame(x1,x2,x3)
colnames(k) <- c("Maximum_Drawdown", "Sharpe_Ratio","Cost_function")
k
p1<-ggplot(data=k,aes(x=Maximum_Drawdown,y=Sharpe_Ratio,color=factor(Cost_function)))+geom_point()+ stat_smooth()
p1
k



plot(0:36,port.simM3(5,0,1000)[,1],type='l',ylim=c(0.9,1.4))
for(i in 1:9){
  lines(0:36,port.simM3(5,0,1000)[,i+1],col=i+1)
}


plot(0:36,port.simM25(25,1000)[,1],type='l',ylim=c(0.9,1.4))
for(i in 1:9){
  lines(0:36,port.simM25(25,1000)[,i+1],col=i+1)
}


> plot(density(rghyp(100000,best[[2]])),xlim=c(-0.5,0.5))
> plot(density(rnorm(100000,0,sqrt(0.001383374))),xlim=c(-0.5,0.5))
density(rghyp(100000,best[[2]])
        x=c(rghyp(100000,best[[2]]),rnorm(100000,0,sqrt(0.001383374)))
        distribution=c(rep('Normal Inverse Gaussian Distribution',100000),rep('Normal Distribution',100000))
        x=data.frame(x,distribution)
pd<-ggplot(data=x,aes(x=x,color=distribution))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(-0.5, 0.5))
pd
