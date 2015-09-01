########################################################################################################################
OnlyGamma=function(){
  p.1000.0<-port.sim2(100,0,500)
  p.100.0<-port.sim2(10,0,500)
  p.20.0<-port.sim2(2,0,500)
  p.10.0<-port.sim2(1,0,500)
  p.5.0<-port.sim2(0.5,0,500)
  p.1.0<-port.sim2(0.1,0,500)
  p.02.0<-port.sim2(0.02,0,500)
  p.01.0<-port.sim2(0.01,0,500)
  p.005.0<-port.sim2(0.005,0,500)
  p.004.0<-port.sim2(0.004,0,500)
  p.003.0<-port.sim2(0.003,0,500)
  p.002.0<-port.sim2(0.002,0,500)
  p.001.0<-port.sim2(0.001,0,500)
  p.0005.0<-port.sim2(0.0005,0,500)
  p.0004.0<-port.sim2(0.0004,0,500)
  p.0001.0<-port.sim2(0.0001,0,500)
  p.00009.0<-port.sim2(9*0.1^5,0,500)
  p.00005.0<-port.sim2(5*0.1^5,0,500)
  p.00001.0<-port.sim2(0.1^5,0,500)
  
  
  
  f.1000.0<-perfYY(p.1000.0,500)
  f.100.0<-perfYY(p.100.0,500)
  f.20.0<-perfYY(p.20.0,500)
  f.10.0<-perfYY(p.10.0,500)
  f.5.0<-perfYY(p.5.0,500)
  f.1.0<-perfYY(p.1.0,500)
  f.02.0<-perfYY(p.02.0,500)
  f.01.0<-perfYY(p.01.0,500)
  f.005.0<-perfYY(p.005.0,500)
  f.004.0<-perfYY(p.004.0,500)
  f.003.0<-perfYY(p.003.0,500)
  f.002.0<-perfYY(p.002.0,500)
  f.001.0<-perfYY(p.001.0,500)
  f.0005.0<-perfYY(p.0005.0,500)
  f.0004.0<-perfYY(p.0004.0,500)
  f.0001.0<-perfYY(p.0001.0,500)
  f.00009.0<-perfYY(p.00009.0,500)
  f.00005.0<-perfYY(p.00005.0,500)
  f.00001.0<-perfYY(p.00001.0,500)
}
r00=rbind(f.1000.0,f.100.0,f.10.0,f.5.0,f.1.0,f.05.0,f.02.0,f.01.0,f.005.0,f.004.0,f.003.0,f.002.0,
          f.001.0,f.0005.0,f.0004.0,f.0001.0,f.00009.0,f.00005.0,f.00001.0)



Onlygammafigure=function(){
  plot(density(p.100.0[301,]),xlim=c(10000,30000))
  lines(density(p.10.0[301,]),col=2)
  lines(density(p.1.0[301,]),col=4)
  lines(density(p.02.0[301,]),col=5)
  lines(density(p.01.0[301,]),col=5)
  lines(density(p.005.0[301,]),col=6)
  lines(density(p.004.0[301,]),col=7)
  lines(density(p.003.0[301,]),col=8)
  lines(density(p.002.0[301,]),col=9)
  lines(density(p.001.0[301,]),col=10)
  lines(density(p.0005.0[301,]),col=11)
  lines(density(p.0001.0[301,]),col=12)
  lines(density(p.00001.0[301,]),col=13)
  
  
  
  x1<-c(p.100.0[301,],p.1.0[301,],p.02.0[301,],p.01.0[301,],p.005.0[301,],
        p.002.0[301,],p.001.0[301,],p.0005.0[301,],p.0001.0[301,],p.00001.0[301,])
  x2=c(rep("10",500),rep("0.1",500),rep("0.02",500),rep("0.01",500),rep("0.005",500),
       rep("0.002",500),rep("0.001",500),rep("0.0005",500),rep("0.0001",500),rep("0.00001",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Gamma")
  p<-ggplot(data=k,aes(x=Dollar,color=Gamma))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Gamma", summarise, grp.mean=mean(Dollar))
  head(mu)
  p<-p+geom_vline(data=mu, aes(xintercept=grp.mean, color=Gamma),linetype="dashed")
  p
}
########################################################################################################################

G0.002lambda=function(){
  p.002.100<-port.sim2(0.002,10,500)
  p.002.10<-port.sim2(0.002,1,500)
  p.002.1<-port.sim2(0.002,0.1,500)
  p.002.2<-port.sim2(0.002,0.1^2,500)
  p.002.3<-port.sim2(0.002,0.1^3,500)
  p.002.4<-port.sim2(0.002,0.1^4,500)
  p.002.5<-port.sim2(0.002,0.1^5,500)
  p.002.6<-port.sim2(0.002,0.1^6,500)
  p.002.7<-port.sim2(0.002,0.1^7,500)
  p.002.8<-port.sim2(0.002,0.1^8,500)
  p.002.9<-port.sim2(0.002,0.1^9,500)
  p.002.10<-port.sim2(0.002,0.1^10,500)
  p.002.11<-port.sim2(0.002,0.1^11,500)
  p.002.13<-port.sim2(0.002,0.1^13,500)
  p.002.15<-port.sim2(0.002,0.1^15,500)
  p.002.18<-port.sim2(0.002,0.1^18,500)
  p.002.23<-port.sim2(0.002,0.1^23,500)
  
  f.002.100<-perfYY(p.002.100,500)
  f.002.10<-perfYY(p.002.10,500)
  f.002.1<-perfYY(p.002.1,500)
  f.002.2<-perfYY(p.002.2,500)
  f.002.3<-perfYY(p.002.3,500)
  f.002.4<-perfYY(p.002.4,500)
  f.002.5<-perfYY(p.002.5,500)
  f.002.6<-perfYY(p.002.6,500)
  f.002.7<-perfYY(p.002.7,500)
  f.002.8<-perfYY(p.002.8,500)
  f.002.9<-perfYY(p.002.9,500)
  f.002.10<-perfYY(p.002.10,500)
  f.002.11<-perfYY(p.002.11,500)
  f.002.12<-perfYY(p.002.12,500)
  f.002.13<-perfYY(p.002.13,500)
  f.002.15<-perfYY(p.002.15,500)
  f.002.18<-perfYY(p.002.18,500)
  f.002.23<-perfYY(p.002.23,500)
}

r02=rbind(f.002.100,f.002.1,f.002.2,f.002.3,f.002.4,f.002.5,f.002.6,f.002.7,f.002.8,f.002.9,f.002.10,f.002.11
          ,f.002.13,f.002.15,f.002.18,f.002.23,f.002.0)

figure<-function(){
  x1<-c(p.100.0[301,],p.002.5[301,],p.002.6[301,],p.002.7[301,],p.002.8[301,],
        p.002.9[301,],p.002.10[301,],p.002.11[301,],p.002.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p3<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 18000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p3<-p3+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p3
}
########################################################################################################################

G0.5lambda=function(){
  p.5.100<-port.sim2(0.5,10,500)
  p.5.10<-port.sim2(0.5,1,500)
  p.5.1<-port.sim2(0.5,0.1,500)
  p.5.2<-port.sim2(0.5,0.1^2,500)
  p.5.3<-port.sim2(0.5,0.1^3,500)
  p.5.4<-port.sim2(0.5,0.1^4,500)
  p.5.5<-port.sim2(0.5,0.1^5,500)
  p.5.6<-port.sim2(0.5,0.1^6,500)
  p.5.7<-port.sim2(0.5,0.1^7,500)
  p.5.8<-port.sim2(0.5,0.1^8,500)
  p.5.9<-port.sim2(0.5,0.1^9,500)
  p.5.10<-port.sim2(0.5,0.1^10,500)
  p.5.11<-port.sim2(0.5,0.1^11,500)
  p.5.13<-port.sim2(0.5,0.1^13,500)
  p.5.15<-port.sim2(0.5,0.1^15,500)
  p.5.18<-port.sim2(0.5,0.1^18,500)
  p.5.23<-port.sim2(0.5,0.1^23,500)
  
  f.5.100<-perfYY(p.5.100,500)
  f.5.10<-perfYY(p.5.10,500)
  f.5.1<-perfYY(p.5.1,500)
  f.5.2<-perfYY(p.5.2,500)
  f.5.3<-perfYY(p.5.3,500)
  f.5.4<-perfYY(p.5.4,500)
  f.5.5<-perfYY(p.5.5,500)
  f.5.6<-perfYY(p.5.6,500)
  f.5.7<-perfYY(p.5.7,500)
  f.5.8<-perfYY(p.5.8,500)
  f.5.9<-perfYY(p.5.9,500)
  f.5.10<-perfYY(p.5.10,500)
  f.5.11<-perfYY(p.5.11,500)
  f.5.12<-perfYY(p.5.12,500)
  f.5.13<-perfYY(p.5.13,500)
  f.5.15<-perfYY(p.5.15,500)
  f.5.18<-perfYY(p.5.18,500)
  f.5.23<-perfYY(p.5.23,500)
}

r05=rbind(f.5.100,f.5.1,f.5.2,f.5.3,f.5.4,f.5.5,f.5.6,f.5.7,f.5.8,f.5.9,f.5.10,f.5.11
          ,f.5.13,f.5.15,f.5.18,f.5.23,f.5.0)

figure<-function(){
  x1<-c(p.100.0[301,],p.5.5[301,],p.5.6[301,],p.5.7[301,],p.5.8[301,],
        p.5.9[301,],p.5.10[301,],p.5.11[301,],p.5.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p1<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 14000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p1<-p1+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p1
}
########################################################################################################################

G0.0001lambda=function(){
  p.0001.100<-port.sim2(0.0001,10,500)
  p.0001.10<-port.sim2(0.0001,1,500)
  p.0001.1<-port.sim2(0.0001,0.1,500)
  p.0001.2<-port.sim2(0.0001,0.1^2,500)
  p.0001.3<-port.sim2(0.0001,0.1^3,500)
  p.0001.4<-port.sim2(0.0001,0.1^4,500)
  p.0001.5<-port.sim2(0.0001,0.1^5,500)
  p.0001.6<-port.sim2(0.0001,0.1^6,500)
  p.0001.7<-port.sim2(0.0001,0.1^7,500)
  p.0001.8<-port.sim2(0.0001,0.1^8,500)
  p.0001.9<-port.sim2(0.0001,0.1^9,500)
  p.0001.10<-port.sim2(0.0001,0.1^10,500)
  p.0001.11<-port.sim2(0.0001,0.1^11,500)
  p.0001.13<-port.sim2(0.0001,0.1^13,500)
  p.0001.15<-port.sim2(0.0001,0.1^15,500)
  p.0001.18<-port.sim2(0.0001,0.1^18,500)
  p.0001.23<-port.sim2(0.0001,0.1^23,500)
  
  f.0001.100<-perfYY(p.0001.100,500)
  f.0001.10<-perfYY(p.0001.10,500)
  f.0001.1<-perfYY(p.0001.1,500)
  f.0001.2<-perfYY(p.0001.2,500)
  f.0001.3<-perfYY(p.0001.3,500)
  f.0001.4<-perfYY(p.0001.4,500)
  f.0001.5<-perfYY(p.0001.5,500)
  f.0001.6<-perfYY(p.0001.6,500)
  f.0001.7<-perfYY(p.0001.7,500)
  f.0001.8<-perfYY(p.0001.8,500)
  f.0001.9<-perfYY(p.0001.9,500)
  f.0001.10<-perfYY(p.0001.10,500)
  f.0001.11<-perfYY(p.0001.11,500)
  f.0001.12<-perfYY(p.0001.12,500)
  f.0001.13<-perfYY(p.0001.13,500)
  f.0001.15<-perfYY(p.0001.15,500)
  f.0001.18<-perfYY(p.0001.18,500)
  f.0001.23<-perfYY(p.0001.23,500)
}
r01=rbind(f.0001.100,f.0001.1,f.0001.2,f.0001.3,f.0001.4,f.0001.5,f.0001.6,f.0001.7,
          f.0001.8,f.0001.9,f.0001.10,f.0001.11,f.0001.13,f.0001.15,f.0001.18,f.0001.23,f.0001.0)


figure0.0001<-function(){
  plot(density(p.100.0[301,]),xlim=c(10000,30000))
  lines(density(p.0001.5[301,]),col=2)
  lines(density(p.0001.6[301,]),col=3)
  lines(density(p.0001.7[301,]),col=4)
  lines(density(p.0001.8[301,]),col=3)
  lines(density(p.0001.9[301,]),col=3)
  lines(density(p.0001.10[301,]),col=5)
  lines(density(p.0001.11[301,]),col=6)
  lines(density(p.0001.13[301,]),col=4)
  lines(density(p.0001.0[301,]),col=2)
  
  x1<-c(p.100.0[301,],p.0001.5[301,],p.0001.6[301,],p.0001.7[301,],p.0001.8[301,],
        p.0001.9[301,],p.0001.10[301,],p.0001.11[301,],p.0001.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p2<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p2<-p2+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p2
}

########################################################################################################################

G0.0005lambda=function(){
  p.0005.100<-port.sim2(0.0005,10,500)
  p.0005.1<-port.sim2(0.0005,0.1,500)
  p.0005.4<-port.sim2(0.0005,0.1^4,500)
  p.0005.5<-port.sim2(0.0005,0.1^5,500)
  p.0005.6<-port.sim2(0.0005,0.1^6,500)
  p.0005.7<-port.sim2(0.0005,0.1^7,500)
  p.0005.85<-port.sim2(0.0005,5*0.1^8,500)
  p.0005.8<-port.sim2(0.0005,0.1^8,500)
  p.0005.95<-port.sim2(0.0005,5*0.1^9,500)
  p.0005.9<-port.sim2(0.0005,0.1^9,500)
  p.0005.10<-port.sim2(0.0005,0.1^10,500)
  p.0005.11<-port.sim2(0.0005,0.1^11,500)
  p.0005.13<-port.sim2(0.0005,0.1^13,500)
  p.0005.15<-port.sim2(0.0005,0.1^15,500)

  
 
  f.0005.100<-perfYY(p.0005.100,500)
  f.0005.1<-perfYY(p.0005.1,500)
  f.0005.4<-perfYY(p.0005.4,500)
  f.0005.5<-perfYY(p.0005.5,500)
  f.0005.6<-perfYY(p.0005.6,500)
  f.0005.7<-perfYY(p.0005.7,500)
  f.0005.85<-perfYY(p.0005.85,500)
  f.0005.8<-perfYY(p.0005.8,500)
  f.0005.95<-perfYY(p.0005.95,500)
  f.0005.9<-perfYY(p.0005.9,500)
  f.0005.10<-perfYY(p.0005.10,500)
  f.0005.11<-perfYY(p.0005.11,500)
  f.0005.13<-perfYY(p.0005.13,500)
  f.0005.15<-perfYY(p.0005.15,500)
}
r0005=rbind(f.0005.100,f.0005.1,f.0005.4,f.0005.5,f.0005.6,f.0005.7,f.0005.85,
          f.0005.8,f.0005.95,f.0005.9,f.0005.10,f.0005.11,f.0005.13,f.0005.15,f.0005.0)


figure0.0005<-function(){
  plot(density(p.100.0[301,]),xlim=c(10000,30000))
  lines(density(p.0001.5[301,]),col=2)
  lines(density(p.0001.6[301,]),col=3)
  lines(density(p.0001.7[301,]),col=4)
  lines(density(p.0001.8[301,]),col=3)
  lines(density(p.0001.9[301,]),col=3)
  lines(density(p.0001.10[301,]),col=5)
  lines(density(p.0001.11[301,]),col=6)
  lines(density(p.0001.13[301,]),col=4)
  lines(density(p.0001.0[301,]),col=2)
  
  x1<-c(p.100.0[301,],p.0005.5[301,],p.0005.6[301,],p.0005.7[301,],p.0005.8[301,],
        p.0005.9[301,],p.0005.10[301,],p.0005.11[301,],p.0005.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p0005<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 25000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p0005<-p0005+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p0005
}