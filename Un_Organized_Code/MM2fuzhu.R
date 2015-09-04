M2OnlyGamma=function(){
  m2p.1000.0<-port.simM2(100,0,500)
  m2p.100.0<-port.simM2(10,0,500)
  m2p.20.0<-port.simM2(2,0,0,500)
  m2p.10.0<-port.simM2(1,0,0,500)
  m2p.5.0<-port.simM2(0.5,0,0,500)
  m2p.1.0<-port.simM2(0.1,0,0,500)
  m2p.02.0<-port.simM2(0.02,0,0,500)
  m2p.01.0<-port.simM2(0.01,0,0,500)
  m2p.005.0<-port.simM2(0.005,0,0,500)
  m2p.004.0<-port.simM2(0.004,0,0,500)
  m2p.003.0<-port.simM2(0.003,0,0,500)
  m2p.002.0<-port.simM2(0.002,0,0,500)
  m2p.001.0<-port.simM2(0.001,0,0,500)
  m2p.0005.0<-port.simM2(0.0005,0,0,500)
  m2p.0004.0<-port.simM2(0.0004,0,0,500)
  m2p.0001.0<-port.simM2(0.0001,0,0,500)
  m2p.00009.0<-port.simM2(9*0.1^5,0,0,500)
  m2p.00005.0<-port.simM2(5*0.1^5,0,0,500)
  m2p.00001.0<-port.simM2(0.1^5,0,0,500)
  
  
  
  m2f.1000.0<-perfYY(m2p.1000.0,500)
  m2f.100.0<-perfYY(m2p.100.0,500)
  m2f.20.0<-perfYY(m2p.20.0,500)
  m2f.10.0<-perfYY(m2p.10.0,500)
  m2f.5.0<-perfYY(m2p.5.0,500)
  m2f.1.0<-perfYY(m2p.1.0,500)
  m2f.02.0<-perfYY(m2p.02.0,500)
  m2f.01.0<-perfYY(m2p.01.0,500)
  m2f.005.0<-perfYY(m2p.005.0,500)
  m2f.004.0<-perfYY(m2p.004.0,500)
  m2f.003.0<-perfYY(m2p.003.0,500)
  m2f.002.0<-perfYY(m2p.002.0,500)
  m2f.001.0<-perfYY(m2p.001.0,500)
  m2f.0005.0<-perfYY(m2p.0005.0,500)
  m2f.0004.0<-perfYY(m2p.0004.0,500)
  m2f.0001.0<-perfYY(m2p.0001.0,500)
  m2f.00009.0<-perfYY(m2p.00009.0,500)
  m2f.00005.0<-perfYY(m2p.00005.0,500)
  m2f.00001.0<-perfYY(m2p.00001.0,500)
}
r00=rbind(m2f.1000.0,m2f.100.0,m2f.10.0,m2f.5.0,m2f.1.0,m2f.02.0,m2f.01.0,m2f.005.0,m2f.004.0,m2f.003.0,m2f.002.0,
          m2f.001.0,m2f.0005.0,m2f.0004.0,m2f.0001.0,m2f.00009.0,m2f.00005.0,m2f.00001.0)

r00

Onlygammafigure=function(){
  plot(density(m2p.100.0[301,]),xlim=c(10000,30000))
  lines(density(m2p.10.0[301,]),col=2)
  lines(density(m2p.1.0[301,]),col=4)
  lines(density(m2p.02.0[301,]),col=5)
  lines(density(m2p.01.0[301,]),col=5)
  lines(density(m2p.005.0[301,]),col=6)
  lines(density(m2p.004.0[301,]),col=7)
  lines(density(m2p.003.0[301,]),col=8)
  lines(density(m2p.002.0[301,]),col=9)
  lines(density(m2p.001.0[301,]),col=10)
  lines(density(m2p.0005.0[301,]),col=11)
  lines(density(m2p.0001.0[301,]),col=12)
  lines(density(m2p.00001.0[301,]),col=13)
  
  
  
  x1<-c(m2p.100.0[301,],m2p.1.0[301,],m2p.02.0[301,],m2p.01.0[301,],m2p.005.0[301,],
        m2p.002.0[301,],m2p.001.0[301,],m2p.0005.0[301,],m2p.0001.0[301,],m2p.00001.0[301,])
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
  m2p.002.100<-port.simM2(0.002,10,0,500)
  m2p.002.z10<-port.simM2(0.002,1,0,500)
  m2p.002.1<-port.simM2(0.002,0.1,0,500)
  m2p.002.2<-port.simM2(0.002,0.1^2,0,500)
  m2p.002.3<-port.simM2(0.002,0.1^3,0,500)
  m2p.002.4<-port.simM2(0.002,0.1^4,0,500)
  m2p.002.5<-port.simM2(0.002,0.1^5,0,500)
  m2p.002.6<-port.simM2(0.002,0.1^6,0,500)
  m2p.002.7<-port.simM2(0.002,0.1^7,0,500)
  m2p.002.8<-port.simM2(0.002,0.1^8,0,500)
  m2p.002.9<-port.simM2(0.002,0.1^9,0,500)
  m2p.002.10<-port.simM2(0.002,0.1^10,0,500)
  m2p.002.11<-port.simM2(0.002,0.1^11,0,500)
  m2p.002.13<-port.simM2(0.002,0.1^13,0,500)
  m2p.002.15<-port.simM2(0.002,0.1^15,0,500)
  m2p.002.18<-port.simM2(0.002,0.1^18,0,500)
  m2p.002.23<-port.simM2(0.002,0.1^23,0,500)
  
  m2f.002.100<-perfYY(m2p.002.100,500)
  m2f.002.z10<-perfYY(m2p.002.z10,500)
  m2f.002.1<-perfYY(m2p.002.1,500)
  m2f.002.2<-perfYY(m2p.002.2,500)
  m2f.002.3<-perfYY(m2p.002.3,500)
  m2f.002.4<-perfYY(m2p.002.4,500)
  m2f.002.5<-perfYY(m2p.002.5,500)
  m2f.002.6<-perfYY(m2p.002.6,500)
  m2f.002.7<-perfYY(m2p.002.7,500)
  m2f.002.8<-perfYY(m2p.002.8,500)
  m2f.002.9<-perfYY(m2p.002.9,500)
  m2f.002.10<-perfYY(m2p.002.10,500)
  m2f.002.11<-perfYY(m2p.002.11,500)
  m2f.002.12<-perfYY(m2p.002.12,500)
  m2f.002.13<-perfYY(m2p.002.13,500)
  m2f.002.15<-perfYY(m2p.002.15,500)
  m2f.002.18<-perfYY(m2p.002.18,500)
  m2f.002.23<-perfYY(m2p.002.23,500)
}

r02=rbind(m2f.002.100,m2f.002.1,m2f.002.2,m2f.002.3,m2f.002.4,m2f.002.5,
          m2f.002.6,m2f.002.7,m2f.002.8,m2f.002.9,m2f.002.10,m2f.002.11
          ,m2f.002.13,m2f.002.15,m2f.002.18,m2f.002.23,m2f.002.0)
r02
figure<-function(){
  x1<-c(m2p.100.0[301,],m2p.002.5[301,],m2p.002.6[301,],m2p.002.7[301,],m2p.002.8[301,],
        m2p.002.9[301,],m2p.002.10[301,],m2p.002.11[301,],m2p.002.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p3<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 30000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p3<-p3+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p3
}
########################################################################################################################

G0.5lambda=function(){
  m2p.5.100<-port.simM2(0.5,10,0,500)
  m2p.5.z10<-port.simM2(0.5,1,0,500)
  m2p.5.1<-port.simM2(0.5,0.1,0,500)
  m2p.5.2<-port.simM2(0.5,0.1^2,0,500)
  m2p.5.3<-port.simM2(0.5,0.1^3,0,500)
  m2p.5.4<-port.simM2(0.5,0.1^4,0,500)
  m2p.5.5<-port.simM2(0.5,0.1^5,0,500)
  m2p.5.6<-port.simM2(0.5,0.1^6,0,500)
  m2p.5.7<-port.simM2(0.5,0.1^7,0,500)
  m2p.5.8<-port.simM2(0.5,0.1^8,0,500)
  m2p.5.9<-port.simM2(0.5,0.1^9,0,500)
  m2p.5.10<-port.simM2(0.5,0.1^10,0,500)
  m2p.5.11<-port.simM2(0.5,0.1^11,0,500)
  m2p.5.13<-port.simM2(0.5,0.1^13,0,500)
  m2p.5.15<-port.simM2(0.5,0.1^15,0,500)
  m2p.5.18<-port.simM2(0.5,0.1^18,0,500)
  m2p.5.23<-port.simM2(0.5,0.1^23,0,500)
  
  m2f.5.100<-perfYY(m2p.5.100,500)
  m2f.5.z10<-perfYY(m2p.5.z10,500)
  m2f.5.1<-perfYY(m2p.5.1,500)
  m2f.5.2<-perfYY(m2p.5.2,500)
  m2f.5.3<-perfYY(m2p.5.3,500)
  m2f.5.4<-perfYY(m2p.5.4,500)
  m2f.5.5<-perfYY(m2p.5.5,500)
  m2f.5.6<-perfYY(m2p.5.6,500)
  m2f.5.7<-perfYY(m2p.5.7,500)
  m2f.5.8<-perfYY(m2p.5.8,500)
  m2f.5.9<-perfYY(m2p.5.9,500)
  m2f.5.10<-perfYY(m2p.5.10,500)
  m2f.5.11<-perfYY(m2p.5.11,500)
  m2f.5.12<-perfYY(m2p.5.12,500)
  m2f.5.13<-perfYY(m2p.5.13,500)
  m2f.5.15<-perfYY(m2p.5.15,500)
  m2f.5.18<-perfYY(m2p.5.18,500)
  m2f.5.23<-perfYY(m2p.5.23,500)
}

r05=rbind(m2f.5.100,m2f.5.1,m2f.5.2,m2f.5.3,m2f.5.4,m2f.5.5,m2f.5.6,m2f.5.7,m2f.5.8,m2f.5.9,m2f.5.10,m2f.5.11
          ,m2f.5.13,m2f.5.15,m2f.5.18,m2f.5.23,m2f.5.0)
r05
figure<-function(){
  x1<-c(m2p.100.0[301,],m2p.5.5[301,],m2p.5.6[301,],m2p.5.7[301,],m2p.5.8[301,],
        m2p.5.9[301,],m2p.5.10[301,],m2p.5.11[301,],m2p.5.0[301,])
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
  m2p.0001.100<-port.simM2(0.0001,10,0,500)
  m2p.0001.z10<-port.simM2(0.0001,1,0,500)
  m2p.0001.1<-port.simM2(0.0001,0.1,0,500)
  m2p.0001.2<-port.simM2(0.0001,0.1^2,0,500)
  m2p.0001.3<-port.simM2(0.0001,0.1^3,0,500)
  m2p.0001.4<-port.simM2(0.0001,0.1^4,0,500)
  m2p.0001.5<-port.simM2(0.0001,0.1^5,0,500)
  m2p.0001.6<-port.simM2(0.0001,0.1^6,0,500)
  m2p.0001.7<-port.simM2(0.0001,0.1^7,0,500)
  m2p.0001.8<-port.simM2(0.0001,0.1^8,0,500)
  m2p.0001.95<-port.simM2(0.0001,5*0.1^9,0,500)
  m2p.0001.9<-port.simM2(0.0001,0.1^9,0,500)
  m2p.0001.105<-port.simM2(0.0001,5*0.1^10,0,500)
  m2p.0001.10<-port.simM2(0.0001,0.1^10,0,500)
  m2p.0001.11<-port.simM2(0.0001,0.1^11,0,500)
  m2p.0001.13<-port.simM2(0.0001,0.1^13,0,500)
  m2p.0001.15<-port.simM2(0.0001,0.1^15,0,500)
  m2p.0001.18<-port.simM2(0.0001,0.1^18,0,500)
  m2p.0001.23<-port.simM2(0.0001,0.1^23,0,500)
  
  m2f.0001.100<-perfYY(m2p.0001.100,500)
  m2f.0001.z10<-perfYY(m2p.0001.z10,500)
  m2f.0001.1<-perfYY(m2p.0001.1,500)
  m2f.0001.2<-perfYY(m2p.0001.2,500)
  m2f.0001.3<-perfYY(m2p.0001.3,500)
  m2f.0001.4<-perfYY(m2p.0001.4,500)
  m2f.0001.5<-perfYY(m2p.0001.5,500)
  m2f.0001.6<-perfYY(m2p.0001.6,500)
  m2f.0001.7<-perfYY(m2p.0001.7,500)
  m2f.0001.8<-perfYY(m2p.0001.8,500)
  m2f.0001.95<-perfYY(m2p.0001.95,500)
  m2f.0001.9<-perfYY(m2p.0001.9,500)
  m2f.0001.105<-perfYY(m2p.0001.105,500)
  m2f.0001.10<-perfYY(m2p.0001.10,500)
  m2f.0001.11<-perfYY(m2p.0001.11,500)
  m2f.0001.12<-perfYY(m2p.0001.12,500)
  m2f.0001.13<-perfYY(m2p.0001.13,500)
  m2f.0001.15<-perfYY(m2p.0001.15,500)
  m2f.0001.18<-perfYY(m2p.0001.18,500)
  m2f.0001.23<-perfYY(m2p.0001.23,500)
}
r01=
  rbind(m2f.0001.100,m2f.0001.1,m2f.0001.2,m2f.0001.3,m2f.0001.4,m2f.0001.5,m2f.0001.6,m2f.0001.7,
          m2f.0001.8,m2f.0001.95,m2f.0001.9,m2f.0001.105,m2f.0001.10,m2f.0001.11,m2f.0001.13,m2f.0001.15,
        m2f.0001.18,m2f.0001.23,m2f.0001.0)

r01
r00


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
  
  x1<-c(m2p.100.0[301,],m2p.0001.5[301,],m2p.0001.6[301,],m2p.0001.7[301,],m2p.0001.8[301,],
        m2p.0001.9[301,],m2p.0001.10[301,],m2p.0001.11[301,],m2p.0001.0[301,])
  x2=c(rep("10",500),rep("1e-05",500),rep("1e-06",500),rep("1e-07",500),rep("1e-08",500),
       rep("1e-09",500),rep("1e-10",500),rep("1e-11",500),rep("0",500))
  k<-data.frame(x1,x2)
  colnames(k) <- c("Dollar", "Lambda")
  p2<-ggplot(data=k,aes(x=Dollar,color=Lambda))+stat_density(geom='line',position='identity')+coord_cartesian(xlim = c(10000, 40000))
  #p<-p+labs(title='Density of Final Wealth For Different Gamma')
  library(plyr)
  mu <- ddply(k, "Lambda", summarise, grp.mean=mean(Dollar))
  head(mu)
  p2<-p2+geom_vline(data=mu, aes(xintercept=grp.mean, color=Lambda),linetype="dashed")
  p2
}
#####################################
M2OnlyGammalambda0.1=function(){
  m2p.1000.8<-port.simM2(100,0.1^8,0,500)
  m2p.100.8<-port.simM2(10,0.1^8,0,500)
  m2p.20.8<-port.simM2(2,0.1^8,0,500)
  m2p.10.8<-port.simM2(1,0.1^8,0,500)
  m2p.5.8<-port.simM2(0.5,0.1^8,0,500)
  m2p.1.8<-port.simM2(0.1,0.1^8,0,500)
  m2p.02.8<-port.simM2(0.02,0.1^8,0,500)
  m2p.01.8<-port.simM2(0.01,0.1^8,0,500)
  m2p.005.8<-port.simM2(0.005,0.1^8,0,500)
  m2p.004.8<-port.simM2(0.004,0.1^8,0,500)
  m2p.003.8<-port.simM2(0.003,0.1^8,0,500)
  m2p.002.8<-port.simM2(0.002,0.1^8,0,500)
  m2p.001.8<-port.simM2(0.001,0.1^8,0,500)
  m2p.0005.8<-port.simM2(0.0005,0.1^8,0,500)
  m2p.0004.8<-port.simM2(0.0004,0.1^8,0,500)
  m2p.0001.8<-port.simM2(0.0001,0.1^8,0,500)
  m2p.00009.8<-port.simM2(9*0.1^5,0.1^8,0,500)
  m2p.00005.8<-port.simM2(5*0.1^5,0.1^8,0,500)
  m2p.00001.8<-port.simM2(0.1^5,0.1^8,0,500)
  
  
  
  m2f.1000.8<-perfYY(m2p.1000.8,500)
  m2f.100.8<-perfYY(m2p.100.8,500)
  m2f.20.8<-perfYY(m2p.20.8,500)
  m2f.10.8<-perfYY(m2p.10.8,500)
  m2f.5.8<-perfYY(m2p.5.8,500)
  m2f.1.8<-perfYY(m2p.1.8,500)
  m2f.02.8<-perfYY(m2p.02.8,500)
  m2f.01.8<-perfYY(m2p.01.8,500)
  m2f.005.8<-perfYY(m2p.005.8,500)
  m2f.004.8<-perfYY(m2p.004.8,500)
  m2f.003.8<-perfYY(m2p.003.8,500)
  m2f.002.8<-perfYY(m2p.002.8,500)
  m2f.001.8<-perfYY(m2p.001.8,500)
  m2f.0005.8<-perfYY(m2p.0005.8,500)
  m2f.0004.8<-perfYY(m2p.0004.8,500)
  m2f.0001.8<-perfYY(m2p.0001.8,500)
  m2f.00009.8<-perfYY(m2p.00009.8,500)
  m2f.00005.8<-perfYY(m2p.00005.8,500)
  m2f.00001.8<-perfYY(m2p.00001.8,500)
}


r008=rbind(m2f.1000.8,m2f.100.8,m2f.10.8,m2f.5.8,m2f.1.8,m2f.02.8,m2f.01.8,m2f.005.8,m2f.004.8,m2f.003.8
          ,m2f.002.8,
          m2f.001.8,m2f.0005.8,m2f.0004.8,m2f.0001.8,m2f.00009.8,m2f.00005.8,m2f.00001.8)
r008
#################################################################################################################
G0.002x=function(){
  m2p.001.0.100<-port.simM2(0.001,0,100,500)
  m2p.001.0.z10<-port.simM2(0.001,0,10,500)
  m2p.001.0.z8<-port.simM2(0.001,0,8,500)
  m2p.001.0.z5<-port.simM2(0.001,0,5,500)
  m2p.001.0.z2<-port.simM2(0.001,0,2,500)
  m2p.001.0.z1<-port.simM2(0.001,0,1,500)
  m2p.001.0.9<-port.simM2(0.001,0,0.9,500)
  m2p.001.0.7<-port.simM2(0.001,0,0.7,500)
  m2p.001.0.5<-port.simM2(0.001,0,0.5,500)
  m2p.001.0.3<-port.simM2(0.001,0,0.3,500)
  m2p.001.0.2<-port.simM2(0.001,0,0.2,500)
  m2p.001.0.1<-port.simM2(0.001,0,0.1,500)
  m2p.001.0.01<-port.simM2(0.001,0,0.01,500)
  m2p.001.0.001<-port.simM2(0.001,0,0.001,500)
  m2p.001.0.001<-port.simM2(0.001,0,0.001,500)
  m2p.001.0.f001<-port.simM2(0.001,0,-0.01,500)
  m2p.001.0.f1<-port.simM2(0.001,0,-0.1,500)
  
  
  m2f.001.0.100<-perfYY(m2p.001.0.100,500)
  m2f.001.0.z10<-perfYY(m2p.001.0.z10,500)
  m2f.001.0.z8<-perfYY(m2p.001.0.z8,500)
  m2f.001.0.z5<-perfYY(m2p.001.0.z5,500)
  m2f.001.0.z2<-perfYY(m2p.001.0.z2,500)
  m2f.001.0.z1<-perfYY(m2p.001.0.z1,500)
  m2f.001.0.9<-perfYY(m2p.001.0.9,500)
  m2f.001.0.7<-perfYY(m2p.001.0.7,500)
  m2f.001.0.5<-perfYY(m2p.001.0.5,500)
  m2f.001.0.3<-perfYY(m2p.001.0.3,500)
  m2f.001.0.2<-perfYY(m2p.001.0.2,500)
  m2f.001.0.1<-perfYY(m2p.001.0.1,500)
  m2f.001.0.01<-perfYY(m2p.001.0.01,500)
  m2f.001.0.001<-perfYY(m2p.001.0.001,500)
  m2f.001.0.f01<-perfYY(m2p.001.0.f001,500)
  m2f.001.0.f1<-perfYY(m2p.001.0.f1,500)
}


rd001=rbind(m2f.001.0.100,m2f.001.0.z10,m2f.001.0.z8,m2f.001.0.z5,m2f.001.0.z2,m2f.001.0.z1
            ,m2f.001.0.7,m2f.001.0.5,m2f.001.0.3,m2f.001.0.2,m2f.001.0.1,m2f.001.0.01,m2f.001.0.001,m2f.001.0.f01,m2f.001.0.f1)
rd001
plot(density(m2p.001.0.100[301,]),ylim=c(0,0.0002))
lines(density(m2p.001.0.z10[301,]))
lines(density(m2p.001.0.z8[301,]))
lines(density(m2p.001.0.z5[301,]))
lines(density(m2p.001.0.z2[301,]))
lines(density(m2p.001.0.z1[301,]))
lines(density(m2p.001.0.9[301,]))
lines(density(m2p.001.0.7[301,]))
lines(density(m2p.001.0.5[301,]))
lines(density(m2p.001.0.3[301,]))
lines(density(m2p.001.0.2[301,]))
lines(density(m2p.001.0.1[301,]))
lines(density(m2p.001.0.01[301,]))
lines(density(m2p.001.0.001[301,]))
lines(density(m2p.001.0[301,]),col=2)

plot(r00[,3],r00[,1],type='l')
lines(rd001[,3],rd001[,1])
################################################################################################################
G0.002lambda=function(){
  m2p.005.0.100<-port.simM2(0.005,0,100,500)
  m2p.005.0.z1<-port.simM2(0.005,0,1,500)
  m2p.005.0.7<-port.simM2(0.005,0,0.7,500)
  m2p.005.0.5<-port.simM2(0.005,0,0.5,500)
  m2p.005.0.3<-port.simM2(0.005,0,0.3,500)
  m2p.005.0.1<-port.simM2(0.005,0,0.1,500)
  m2p.005.0.01<-port.simM2(0.005,0,0.01,500)
  m2p.005.0.001<-port.simM2(0.005,0,0.001,500)
  
  m2f.005.0.100<-perfYY(m2p.005.0.100,500)
  m2f.005.0.z1<-perfYY(m2p.005.0.z1,500)
  m2f.005.0.7<-perfYY(m2p.005.0.7,500)
  m2f.005.0.5<-perfYY(m2p.005.0.5,500)
  m2f.005.0.3<-perfYY(m2p.005.0.3,500)
  m2f.005.0.1<-perfYY(m2p.005.0.1,500)
  m2f.005.0.01<-perfYY(m2p.005.0.01,500)
  m2f.005.0.001<-perfYY(m2p.005.0.001,500)
}

rd005=rbind(m2f.005.0.100,m2f.005.0.z1,m2f.005.0.7,m2f.005.0.5,m2f.005.0.3,m2f.005.0.1,m2f.005.0.01,m2f.005.0.001,m2f.005.0)

m2p.02.0.100<-port.simM2(0.02,0,100,500)
m2p.02.0.z1<-port.simM2(0.02,0,1,500)
m2p.02.0.7<-port.simM2(0.02,0,0.7,500)
m2p.02.0.5<-port.simM2(0.02,0,0.5,500)
m2p.02.0.3<-port.simM2(0.02,0,0.3,500)
m2p.02.0.1<-port.simM2(0.02,0,0.1,500)
m2p.02.0.01<-port.simM2(0.02,0,0.01,500)
m2p.02.0.001<-port.simM2(0.02,0,0.001,500)

m2f.02.0.100<-perfYY(m2p.02.0.100,500)
m2f.02.0.z1<-perfYY(m2p.02.0.z1,500)
m2f.02.0.7<-perfYY(m2p.02.0.7,500)
m2f.02.0.5<-perfYY(m2p.02.0.5,500)
m2f.02.0.3<-perfYY(m2p.02.0.3,500)
m2f.02.0.1<-perfYY(m2p.02.0.1,500)
m2f.02.0.01<-perfYY(m2p.02.0.01,500)
m2f.02.0.001<-perfYY(m2p.02.0.001,500)
rd02=rbind(m2f.02.0.100,m2f.02.0.z1,m2f.02.0.7,m2f.02.0.5,m2f.02.0.3,m2f.02.0.1,m2f.02.0.01,m2f.02.0.001)


r00

plot(r00[,3],r00[,1],type='l')
lines(sen[,3],sen[,1],col=4)
lines(r02[,3],r02[,1],col=3,lwd=2)
lines(rd02[,3],rd02[,1],col=2)
lines(rd001[,8],rd001[,1],col=5)



sen=rbind(fs2,fs4,fs6,fs9,fs12,fs16,fs20,fs25,fs30,fs40,fs50,fs70)
g0=rbind(ft0,ft1,ft2,ft4,ft6,ft9,ft10,ft12,ft15,ft16,ft20,ft25,ft30,ft40,ft50,ft70,ft90,ft100)
l0=rbind(ft20g1,ft20g2,ft20g4,ft20g8,ft20g20,ft20g30,ft20g50,ft20g80,ft20g400,ft20g1000)
l01=c(0.026,0.026,0.026,0.026,0.025,0.025,0.024,0.024)
l02=c(0.111,0.111,0.111,0.109,0.106,0.095,0.087,0.079)
t20=rbind(ft20,ft20g0t0.1,ft20g0t1,ft20g0t2,ft20g0t4,ft20g0t6,ft20g0t10,ft20g0t15)
t1=c(0.026,0.033,0.035,0.046,0.065,0.083,0.090,0.097)
t2=c(0.111,0.093,0.052,0.022,-0.010,-0.023,-0.019,-0.013)
sen

plot(g0[,3],g0[,1],type='l')
lines(sen[,3],sen[,1],col=4)
lines(l0[,3],l0[,1],col=3,lwd=2)
lines(l01,l02,col=2)
lines(t20[,3],t20[,1],col=5)
lines(t1,t2,col=3)
lines(KKKK[,3],KKKK[,1],col=6)
t20





gt20<-port.simM2(6,20,500)
fgt20<-perfYY(gt20,500)
gt80<-port.simM2(6,80,500)
fgt80<-perfYY(gt80,500)
gt400<-port.simM2(6,400,500)
fgt400<-perfYY(gt400,500)
gt1000<-port.simM2(6,1000,500)
fgt1000<-perfYY(gt1000,500)
gt5000<-port.simM2(6,5000,500)
fgt5000<-perfYY(gt5000,500)
gt10000<-port.simM2(6,10000,500)
fgt10000<-perfYY(gt10000,500)
gt20000<-port.simM2(6,20000,500)
fgt20000<-perfYY(gt20000,500)



KKKK=rbind(ft6,fgt20,fgt80,fgt400,fgt1000,fgt5000,fgt10000,fgt20000)

length(sen[,3])

x1<-as.numeric(c(g0[,3],KKKK[,3],l0[,3],sen[,3]))
x2<-as.numeric(c(g0[,1],KKKK[,1],l0[,1],sen[,1]))
x3=c(rep('Mean-Variance',18),rep('MVK lambda 6',8),rep('MVK lambda 20',10),rep('Risk sensitive',12))

k<-data.frame(x1,x2,x3)
colnames(k) <- c("Maximum_Drawdown", "Sharpe_Ratio","Cost_function")
k
p1<-ggplot(data=k,aes(x=Maximum_Drawdown,y=Sharpe_Ratio,color=factor(Cost_function)))+geom_point()+ stat_smooth()
p1
KKKK[,c(8,9,1,3,4,5,7)]
