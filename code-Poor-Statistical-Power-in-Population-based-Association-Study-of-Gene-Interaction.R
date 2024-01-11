#library related packages
library(reshape2) 
library(ggplot2)   
library(ggpubr)

#basic function
#Pm1,Pm2,Px1,Px2 are the probabilities of carrying gene M1, M2, X1, and X2 by the affected individual,respectively.
#d1 is the LD between X1 and M1, d2 is the LD between X2 and M2, k is the coefficient of interaction, r1 is the relative risks of X1, and r2 is the relative risks of X2.
#PD is the disease prevalence, nD is the size of N chromosomes for the case-control sample, sig.level is the significant criterion p-value, and df is the degree of freedom.
myfunctionE<-function(Pm1,Pm2,Px1,Px2,d1,d2,k,r1,r2,PD,nD,sig.level,df){   
  #Calculate frequency of chromosomes in case samples that carrying X1, X2,The notation * indicates any allele of 0 or 1 at the current site.fx1d is f1*D, fx2d is f*1D.
  fx1d<-r1*Px1/(1+(r1-1)*Px1)           
  fx2d<-r2*Px2/(1+(r2-1)*Px2)    
  #Let nd equals nD
  nd<-nD
  #Calculate f11D, f10D, f01D, and f00D in case samples, f notes the frequencies of disease genes (X1 and X2)
  f11D<-(1+k)*fx1d*fx2d                
  f10D<-fx1d-(1+k)*fx1d*fx2d            
  f01D<-fx2d-(1+k)*fx1d*fx2d            
  f00D<-1-fx1d-fx2d+(1+k)*fx1d*fx2d  
  #Calculate t1,t2 (t1,t2 are parameters required to calculate frequencies in the control group)
  ta<-1-fx1d*PD/Px1                    
  tb<-1-(1-fx1d)*PD/(1-Px1)
  tc<-1-fx2d*PD/Px2
  td<-1-(1-fx2d)*PD/(1-Px2)
  t1<-ta/tb
  t2<-tc/td
  #Calculate Calculate frequency of chromosomes in control samples that carrying X1, X2, The notation * indicates any allele of 0 or 1 at the current site.dfx1 is f1*d, dfx2 is f*1d.
  dfx1<-t1*Px1/(1+(t1-1)*Px1)          
  dfx2<-t2*Px2/(1+(t2-1)*Px2)    
  #Calculate the frequencies of disease genes in control samples
  df11<-(1+k)*dfx1*dfx2                 
  df10<-dfx1-(1+k)*dfx1*dfx2            
  df01<-dfx2-(1+k)*dfx1*dfx2           
  df00<-1-dfx1-dfx2+(1+k)*dfx1*dfx2  
  #Calculate P(m1=1|x1=1),P(m2=1|x2=1),P(m1=1|x1=0),P(m2=1|x2=0),P(m1=0|x1=1),P(m2=0|x2=1),P(m1=0|x1=0), and P(m2=0|x2=0)
  Pa1<-(Pm1*Px1+d1)/Px1                
  Pa2<-(Pm2*Px2+d2)/Px2                 
  Pb1<-(Pm1*(1-Px1)+d1)/(1-Px1)         
  Pb2<-(Pm2*(1-Px2)+d2)/(1-Px2)         
  Pc1<-((1-Pm1)*Px1+d1)/Px1             
  Pc2<-((1-Pm2)*Px2+d2)/Px2            
  Pd1<-((1-Pm1)*(1-Px1)+d1)/(1-Px1)    
  Pd2<-((1-Pm2)*(1-Px2)+d2)/(1-Px2)    
  #Calculate 4 frequencies p11D,p22D,p12D, and p21D for the case samples
  p11D<-Pa1*Pa2*f11D+Pb1*Pa2*f01D+Pa1*Pb2*f10D+Pb1*Pb2*f00D
  p22D<-Pc1*Pc2*f11D+Pd1*Pc2*f01D+Pc1*Pd2*f10D+Pd1*Pd2*f00D
  p12D<-Pa1*Pc2*f11D+Pb1*Pc2*f01D+Pa1*Pd2*f10D+Pb1*Pd2*f00D
  p21D<-Pc1*Pa2*f11D+Pd1*Pa2*f01D+Pc1*Pb2*f10D+Pd1*Pb2*f00D
  #Calculate 4 frequencies p11d,p22d,p12d, and p21d for the control samples
  p11d<-Pa1*Pa2*df11+Pb1*Pa2*df01+Pa1*Pb2*df10+Pb1*Pb2*df00 
  p22d<-Pc1*Pc2*df11+Pd1*Pc2*df01+Pc1*Pd2*df10+Pd1*Pd2*df00
  p12d<-Pa1*Pc2*df11+Pb1*Pc2*df01+Pa1*Pd2*df10+Pb1*Pd2*df00
  p21d<-Pc1*Pa2*df11+Pd1*Pa2*df01+Pc1*Pb2*df10+Pd1*Pb2*df00
  #Calculate the expectation E1 of the LD between the marker genes M1 and M2 for the case group and E2 for the control group 
  E1<-p11D*p22D-p12D*p21D              
  E2<-p11d*p22d-p12d*p21d    
  #Calculate the expectation of effect size of gene-gene interaction in the case-control study
  dE<-E1-E2 
  #Calculate the covariance V1 of the difference of frequencies for the case group and the covariance V2 of the difference of frequencies for the control group
  V1<-((p11D+p12D)*(1-p11D-p12D)*(p11D+p21D)*(1-p11D-p21D)+(1-2*p11D-2*p12D)*(1-2*p11D-2*p21D)*(p11D*p22D-p12D*p21D)-(p11D*p22D-p12D*p21D)*(p11D*p22D-p12D*p21D))/nD 
  V2<-((p11d+p12d)*(1-p11d-p12d)*(p11d+p21d)*(1-p11d-p21d)+(1-2*p11d-2*p12d)*(1-2*p11d-2*p21d)*(p11d*p22d-p12d*p21d)-(p11d*p22d-p12d*p21d)*(p11d*p22d-p12d*p21d))/nd     
  #Calculate V, the sum of V1 and V2
  V<-V1+V2  
  #Calculate the statistic T
  T<-dE*dE/V
  #Connect sig.level, df, and power1
  k1<- qchisq(sig.level, df = df, lower = FALSE)
  #Calculate the power1 value in the chi-square test
  power1<-pchisq(k1, df = df, ncp = T, lower = FALSE)
  #Output dE,V,and power1
  result <- data.frame(dE=dE,V=V,power1=power1)
  print(result)
}

myfunctionE2<-function(Pm1,Pm2,Px1,Px2,d1a,d2a,k,r1,r2,PD,nD,nd,sig.level,df){   #d1a=D1*D1,d2a=D2*D2
  #Calculate d1max and d2max
  d1max<-min(Px1*(1-Pm1),Pm1*(1-Px1))           
  d2max<-min(Px2*(1-Pm2),Pm2*(1-Px2))
  #Calculate Lewontin's d1 and Lewontin's d2
  d1<-d1max*d1a          
  d2<-d2max*d2a           
  #Calculate frequency of chromosomes in case samples that carrying X1, X2,The notation * indicates any allele of 0 or 1 at the current site.fx1d is f1*D, fx2d is f*1D.
  fx1d<-r1*Px1/(1+(r1-1)*Px1)           
  fx2d<-r2*Px2/(1+(r2-1)*Px2)           
  #Calculate f11D, f10D, f01D, and f00D in case samples, f notes the frequencies of disease genes (X1 and X2)
  f11D<-(1+k)*fx1d*fx2d                
  f10D<-fx1d-(1+k)*fx1d*fx2d            
  f01D<-fx2d-(1+k)*fx1d*fx2d            
  f00D<-1-fx1d-fx2d+(1+k)*fx1d*fx2d  
  #Calculate t1,t2 (t1,t2 are parameters required to calculate frequencies in the control group)
  ta<-1-fx1d*PD/Px1                    
  tb<-1-(1-fx1d)*PD/(1-Px1)
  tc<-1-fx2d*PD/Px2
  td<-1-(1-fx2d)*PD/(1-Px2)
  t1<-ta/tb
  t2<-tc/td
  #Calculate Calculate frequency of chromosomes in control samples that carrying X1, X2, The notation * indicates any allele of 0 or 1 at the current site.dfx1 is f1*d, dfx2 is f*1d.
  dfx1<-t1*Px1/(1+(t1-1)*Px1)          
  dfx2<-t2*Px2/(1+(t2-1)*Px2)    
  #Calculate the frequencies of disease genes in control samples
  df11<-(1+k)*dfx1*dfx2                 
  df10<-dfx1-(1+k)*dfx1*dfx2            
  df01<-dfx2-(1+k)*dfx1*dfx2           
  df00<-1-dfx1-dfx2+(1+k)*dfx1*dfx2  
  #Calculate P(m1=1|x1=1),P(m2=1|x2=1),P(m1=1|x1=0),P(m2=1|x2=0),P(m1=0|x1=1),P(m2=0|x2=1),P(m1=0|x1=0), and P(m2=0|x2=0)
  Pa1<-(Pm1*Px1+d1)/Px1                
  Pa2<-(Pm2*Px2+d2)/Px2                 
  Pb1<-(Pm1*(1-Px1)+d1)/(1-Px1)         
  Pb2<-(Pm2*(1-Px2)+d2)/(1-Px2)         
  Pc1<-((1-Pm1)*Px1+d1)/Px1             
  Pc2<-((1-Pm2)*Px2+d2)/Px2            
  Pd1<-((1-Pm1)*(1-Px1)+d1)/(1-Px1)    
  Pd2<-((1-Pm2)*(1-Px2)+d2)/(1-Px2)    
  #Calculate 4 frequencies p11D,p22D,p12D, and p21D for the case samples
  p11D<-Pa1*Pa2*f11D+Pb1*Pa2*f01D+Pa1*Pb2*f10D+Pb1*Pb2*f00D
  p22D<-Pc1*Pc2*f11D+Pd1*Pc2*f01D+Pc1*Pd2*f10D+Pd1*Pd2*f00D
  p12D<-Pa1*Pc2*f11D+Pb1*Pc2*f01D+Pa1*Pd2*f10D+Pb1*Pd2*f00D
  p21D<-Pc1*Pa2*f11D+Pd1*Pa2*f01D+Pc1*Pb2*f10D+Pd1*Pb2*f00D
  #Calculate 4 frequencies p11d,p22d,p12d, and p21d for the control samples
  p11d<-Pa1*Pa2*df11+Pb1*Pa2*df01+Pa1*Pb2*df10+Pb1*Pb2*df00 
  p22d<-Pc1*Pc2*df11+Pd1*Pc2*df01+Pc1*Pd2*df10+Pd1*Pd2*df00
  p12d<-Pa1*Pc2*df11+Pb1*Pc2*df01+Pa1*Pd2*df10+Pb1*Pd2*df00
  p21d<-Pc1*Pa2*df11+Pd1*Pa2*df01+Pc1*Pb2*df10+Pd1*Pb2*df00
  #Calculate the expectation E1 of the LD between the marker genes M1 and M2 for the case group and E2 for the control group 
  E1<-p11D*p22D-p12D*p21D              
  E2<-p11d*p22d-p12d*p21d    
  #Calculate the expectation of effect size of gene-gene interaction in the case-control study
  dE<-E1-E2 
  #Calculate the covariance V1 of the difference of frequencies for the case group and the covariance V2 of the difference of frequencies for the control group
  V1<-((p11D+p12D)*(1-p11D-p12D)*(p11D+p21D)*(1-p11D-p21D)+(1-2*p11D-2*p12D)*(1-2*p11D-2*p21D)*(p11D*p22D-p12D*p21D)-(p11D*p22D-p12D*p21D)*(p11D*p22D-p12D*p21D))/nD 
  V2<-((p11d+p12d)*(1-p11d-p12d)*(p11d+p21d)*(1-p11d-p21d)+(1-2*p11d-2*p12d)*(1-2*p11d-2*p21d)*(p11d*p22d-p12d*p21d)-(p11d*p22d-p12d*p21d)*(p11d*p22d-p12d*p21d))/nd     
  #Calculate V, the sum of V1 and V2
  V<-V1+V2  
  #Calculate the statistic T
  T<-dE*dE/V
  #Connect sig.level, df, and powe1r
  k1<- qchisq(sig.level, df = df, lower = FALSE)
  #Calculate the power value in the chi-square test
  power1<-pchisq(k1, df = df, ncp = T, lower = FALSE)
  #Output d1,d2,r1a,r2a,and power
  result <- data.frame(d1=d1,d2=d2,power1=power1)
  print(result)
}

myfunctionE3<-function(Pm1,Pm2,Px1,Px2,d1a,d2a,k,r1,r2,PD,nD,nd,sig.level,df){   #d1a=D1*D1,d2a=D2*D2
  #Calculate d1max and d2max
  d1max<-min(Px1*(1-Pm1),Pm1*(1-Px1))           
  d2max<-min(Px2*(1-Pm2),Pm2*(1-Px2))
  #Calculate Lewontin's d1 and Lewontin's d2
  d1<-d1max*d1a          
  d2<-d2max*d2a           
  #Calculate frequency of chromosomes in case samples that carrying X1, X2,The notation * indicates any allele of 0 or 1 at the current site.fx1d is f1*D, fx2d is f*1D.
  fx1d<-r1*Px1/(1+(r1-1)*Px1)           
  fx2d<-r2*Px2/(1+(r2-1)*Px2)           
  #Calculate f11D, f10D, f01D, and f00D in case samples, f notes the frequencies of disease genes (X1 and X2)
  f11D<-(1+k)*fx1d*fx2d                
  f10D<-fx1d-(1+k)*fx1d*fx2d            
  f01D<-fx2d-(1+k)*fx1d*fx2d            
  f00D<-1-fx1d-fx2d+(1+k)*fx1d*fx2d  
  #Calculate t1,t2 (t1,t2 are parameters required to calculate frequencies in the control group)
  ta<-1-fx1d*PD/Px1                    
  tb<-1-(1-fx1d)*PD/(1-Px1)
  tc<-1-fx2d*PD/Px2
  td<-1-(1-fx2d)*PD/(1-Px2)
  t1<-ta/tb
  t2<-tc/td
  #Calculate Calculate frequency of chromosomes in control samples that carrying X1, X2, The notation * indicates any allele of 0 or 1 at the current site.dfx1 is f1*d, dfx2 is f*1d.
  dfx1<-t1*Px1/(1+(t1-1)*Px1)          
  dfx2<-t2*Px2/(1+(t2-1)*Px2)    
  #Calculate the frequencies of disease genes in control samples
  df11<-(1+k)*dfx1*dfx2                 
  df10<-dfx1-(1+k)*dfx1*dfx2            
  df01<-dfx2-(1+k)*dfx1*dfx2           
  df00<-1-dfx1-dfx2+(1+k)*dfx1*dfx2  
  #Calculate P(m1=1|x1=1),P(m2=1|x2=1),P(m1=1|x1=0),P(m2=1|x2=0),P(m1=0|x1=1),P(m2=0|x2=1),P(m1=0|x1=0), and P(m2=0|x2=0)
  Pa1<-(Pm1*Px1+d1)/Px1                
  Pa2<-(Pm2*Px2+d2)/Px2                 
  Pb1<-(Pm1*(1-Px1)+d1)/(1-Px1)         
  Pb2<-(Pm2*(1-Px2)+d2)/(1-Px2)         
  Pc1<-((1-Pm1)*Px1+d1)/Px1             
  Pc2<-((1-Pm2)*Px2+d2)/Px2            
  Pd1<-((1-Pm1)*(1-Px1)+d1)/(1-Px1)    
  Pd2<-((1-Pm2)*(1-Px2)+d2)/(1-Px2)    
  #Calculate 4 frequencies p11D,p22D,p12D, and p21D for the case samples
  p11D<-Pa1*Pa2*f11D+Pb1*Pa2*f01D+Pa1*Pb2*f10D+Pb1*Pb2*f00D
  p22D<-Pc1*Pc2*f11D+Pd1*Pc2*f01D+Pc1*Pd2*f10D+Pd1*Pd2*f00D
  p12D<-Pa1*Pc2*f11D+Pb1*Pc2*f01D+Pa1*Pd2*f10D+Pb1*Pd2*f00D
  p21D<-Pc1*Pa2*f11D+Pd1*Pa2*f01D+Pc1*Pb2*f10D+Pd1*Pb2*f00D
  #Calculate 4 frequencies p11d,p22d,p12d, and p21d for the control samples
  p11d<-Pa1*Pa2*df11+Pb1*Pa2*df01+Pa1*Pb2*df10+Pb1*Pb2*df00 
  p22d<-Pc1*Pc2*df11+Pd1*Pc2*df01+Pc1*Pd2*df10+Pd1*Pd2*df00
  p12d<-Pa1*Pc2*df11+Pb1*Pc2*df01+Pa1*Pd2*df10+Pb1*Pd2*df00
  p21d<-Pc1*Pa2*df11+Pd1*Pa2*df01+Pc1*Pb2*df10+Pd1*Pb2*df00
  #Calculate the expectation E3 of the odd ratio for the case group and the control group 
  E3<-(p11D+p12D)*(p21d+p22d)/((p21D+p22D)*(p11d+p12d))    
  #Calculate the variance of the estimated variance of the logarithm of the odd ratio
  V3<-1/(nD*(p11D+p12D))+1/(nD*(p21D+p22D))+1/(nd*(p11d+p12d))+1/(nd*(p21d+p22d))
  #Calculate the statistic T2
  T2<-log(E3)*log(E3)/V3
  #Connect sig.level, df, and power2
  k1<- qchisq(sig.level, df = df, lower = FALSE)
  #Calculate the power value in the chi-square test
  power2<-pchisq(k1, df = df, ncp = T2, lower = FALSE)
  #Output d1,d2,r1a,r2a,and power2
  result <- data.frame(d1=d1,d2=d2,power2=power2)
  print(result)
}


##Generate data and draw the figure of fig.1
#investigate the relationship between dE and k under different values of r1,r2(r1=r2)
output1a <- data.frame(K=seq(0,5,0.25),
                       Rx=rep(c(2.5,2,1.5),each=21))
tmp <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,output1a$K,output1a$Rx,output1a$Rx,0.1,10000,0.05,1)
output1a$Size <- tmp$dE

data <- output1a
data$Rx<-as.factor(data$Rx)
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=K,y=Size,group=Rx,color=Rx),linewidth=1.3)+
  geom_point(data = data2,aes(x=K,y=Size,group=Rx,color=Rx,shape=Rx),size=3)+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("K")+
  ylab("Effect Size")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0,5.2),breaks = seq(0,5,0.5))+
  scale_y_continuous(limits = c(0,0.04),breaks = seq(0,0.04,0.01))

ggsave("size-K.pdf", plot = p1, height = 4, width = 5)

#investigate the relationship between mathematical expectation (dE) , the relative risk ratio of X1 (r1,r2) under different values of the gene interaction strength of X1 and X2 (k)
output1b <- data.frame(K=rep(c(2,1,0.5),each=21),
                       Rx=seq(1,4,0.15))
tmp <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,output1b$K,output1b$Rx,output1b$Rx,0.1,10000,0.05,1)
output1b$Size <- tmp$dE

data <- output1b
data$K<-as.factor(data$K)
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=Rx,y=Size,group=K,color=K),linewidth=1.3)+
  geom_point(data = data2,aes(x=Rx,y=Size,group=K,color=K,shape=K),size=3)+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("R")+
  ylab("Effect Size")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(1,4.1),breaks = seq(1,4,0.25))+
  scale_y_continuous(limits = c(0,0.04),breaks = seq(0,0.04,0.01))

ggsave("size-RR.pdf", plot = p1, height = 4, width = 5)


##Generate data and draw the figure of fig.2
#investigate the relationship between variance and the total chromosome number of the disease population (nD) under different values of the frequency of allele X1 (Px1)  
output1c <- data.frame(n=seq(3450,50000,2450),
                       Px1=rep(c(0.05,0.1,0.15),each=20))
tmp <- myfunctionE(0.05,0.05,output1c$Px1,0.05,0.0475,0.0475,1,2,2,0.1,output1c$n,0.05,1)
output1c$Variance <- tmp$V

data <- output1c
data$Px1<-as.factor(data$Px1)
data2<-data[-seq(2,59,2),]

p1<-ggplot()+
  geom_line(data = data,aes(x=n,y=Variance,group=Px1,color=Px1),linewidth=0.5)+
  geom_point(data = data2,aes(x=n,y=Variance,group=Px1,color=Px1,shape=Px1),size=2)+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("N")+
  ylab("variance")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0,51000),breaks = seq(0,50000,5000))

ggsave("variance-n.pdf", plot = p1, height = 3, width = 4)


##Generate data and draw the figure of fig.3
#investigate the relationship between power1 and r1,r2(r1=r2) under different values of k(fig.3A)
output1d <- data.frame(Rx=seq(1,4,0.15),
                       k=rep(c(2,1,0.5),each=21))
tmp <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,output1d$k,output1d$Rx,output1d$Rx,0.1,10000,0.05,1)
output1d$power <- tmp$power1

data <- output1d
data$k<-factor(data$k,levels=c("2","1","0.5"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=Rx,y=power,group=k,color=k),linewidth=1.3)+
  geom_point(data = data2,aes(x=Rx,y=power,group=k,color=k,shape=k),size=3)+ 
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+
  xlab("R")+
  ylab("Power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0.9,4.1),breaks = seq(1,4.1,0.2))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))
ggsave("power-r-k.pdf", plot = p1, height = 4, width = 5)

#investigate the relationship between power1 and r1,r2(r1=r2) under different values of nD(fig.3B)
output1e <- data.frame(Rx=seq(1,4,0.15),
                       n=rep(c(10000,5000,1000),each=21))
tmp <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,output1e$Rx,output1e$Rx,0.1,output1e$n,0.05,1)
output1e$power <- tmp$power1

data <- output1e
data$n<-factor(data$n,levels=c("10000","5000","1000"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=Rx,y=power,group=n,color=n),linewidth=1.3)+
  geom_point(data = data2,aes(x=Rx,y=power,group=n,color=n,shape=n),size=3)+ 
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+
  xlab("n")+
  ylab("Power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0.9,4.1),breaks = seq(1,4.1,0.2))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))

ggsave("power-r-n.pdf", plot = p1, height = 4, width = 5)


##Generate data and draw the figure of fig.4
#investigate the relationship between power1 and k under different values of r1,r2(r1=r2)(fig.4A)
output1f <- data.frame(K=seq(0,5,0.25),
                       Rx=rep(c(2.5,2,1.2),each=21))
tmp <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,output1f$K,output1f$Rx,output1f$Rx,0.1,10000,0.05,1)
output1f$power <- tmp$power1

data <- output1f
data$Rx<-factor(data$Rx,levels=c("2.5","2","1.2"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=K,y=power,group=Rx,color=Rx),linewidth=1.3)+
  geom_point(data = data2,aes(x=K,y=power,group=Rx,color=Rx,shape=Rx),size=3)+ 
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("K")+
  ylab("Power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0,5),breaks = seq(0,5,0.5))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))
ggsave("power-k-r.pdf", plot = p1, height = 4, width = 5)

#investigate the relationship between power1 and k under different values of nD(fig.4B)
output1g <- data.frame(K=seq(0,5,0.25),
                       n=rep(c(10000,5000,1000),each=21))
tmp <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,output1g$K,2,2,0.1,output1g$n,0.05,1)
output1g$power <- tmp$power1

data <- output1g
data$n<-factor(data$n,levels=c("10000","5000","1000"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=K,y=power,group=n,color=n),linewidth=1.3)+
  geom_point(data = data2,aes(x=K,y=power,group=n,color=n,shape=n),size=3)+ 
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("K")+
  ylab("Power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0,5),breaks = seq(0,5,0.5))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))
ggsave("power-K-n.pdf", plot = p1, height = 4, width = 5)


##Generate data and draw the figure of fig.5
#investigate the relationship between power1,power2 and d1a (d1/d1max) under different values of Px1 and level of significance (sig.level)
#Generate data of fig.5
output1h <- data.frame(Px1=rep(c(0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.1),each=20),
                       d1a=seq(0.05,1,0.05),
                       siglevel=rep(c(0.05,0.0005,0.000005,0.00000005),each=60))

c1<-list()
for (i in 1:length(output1h$Px1)) {
  tmp1<-myfunctionE2(0.2,0.2,output1h$Px1[i],0.2,output1h$d1a[i],1,1,2,2,0.1,10000,10000,output1h$siglevel[i],1)
  pow<-tmp1$power1
  c1[i]<-pow
}
power1<-as.data.frame(do.call(rbind, c1))
output1h$power1 <- power1$V1

d1<-list()
for (i in 1:length(output1h$Px1)) {
  tmp1<-myfunctionE3(0.2,0.2,output1h$Px1[i],0.2,output1h$d1a[i],1,1,2,2,0.1,10000,10000,output1h$siglevel[i],1)
  pow<-tmp1$power2
  d1[i]<-pow
}
power2<-as.data.frame(do.call(rbind, d1))
output1h$power2 <- power2$V1

#draw fig5A
data <- output1h[1:60,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,1.5),seq(22,38,1.5),seq(42,58,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("d1a", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("d1a", "Px1"),
                                     names_to = "group", values_to = "power") 

p1<-ggplot()+
  geom_line(data=data_long,aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=1.3)+
  geom_point(data=data_long2,aes(x = d1a, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("d1a")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#draw fig5B
data <- output1h[61:120,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,1.5),seq(22,38,1.5),seq(42,58,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("d1a", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("d1a", "Px1"),
                                     names_to = "group", values_to = "power") 
p2<-ggplot()+
  geom_line(data=data_long,aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=1.3)+
  geom_point(data=data_long2,aes(x = d1a, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("d1a")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#draw fig5C
data <- output1h[121:180,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,1.5),seq(22,38,1.5),seq(42,58,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("d1a", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("d1a", "Px1"),
                                     names_to = "group", values_to = "power") 
p3<-ggplot()+
  geom_line(data=data_long,aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=1.3)+
  geom_point(data=data_long2,aes(x = d1a, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("d1a")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#draw fig5D
data <- output1h[181:240,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,1.5),seq(22,38,1.5),seq(42,58,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("d1a", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("d1a", "Px1"),
                                     names_to = "group", values_to = "power") 
p4<-ggplot()+
  geom_line(data=data_long,aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=1.3)+
  geom_point(data=data_long2,aes(x = d1a, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("d1a")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#save fig.5A, fig.5B, fig.5C, and fig.5D
ggsave("power,power2&dla-px1-1.pdf", plot = p1, height = 4, width = 5)
ggsave("power,power2&dla-px1-2.pdf", plot = p2, height = 4, width = 5)
ggsave("power,power2&dla-px1-3.pdf", plot = p3, height = 4, width = 5)
ggsave("power,power2&dla-px1-4.pdf", plot = p4, height = 4, width = 5)


##Generate data and draw the figure of fig.6
#investigate the relationship between power1,power2 and d1a under different values of Px1 and k
#Generate data of fig.6
output1j <- data.frame(Px1=rep(c(0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.1),each=21),
                       k=seq(0,2,0.1),
                       siglevel=rep(c(0.05,0.0005,0.000005,0.00000005),each=63))
c1<-list()
for (i in 1:length(output1j$Px1)) {
  tmp1<-myfunctionE2(0.2,0.2,output1j$Px1[i],0.2,0.5,1,output1j$k[i],2,2,0.1,10000,10000,output1j$siglevel[i],1)
  pow<-tmp1$power1
  c1[i]<-pow
}
power1<-as.data.frame(do.call(rbind, c1))
output1j$power1 <- power1$V1

d1<-list()
for (i in 1:length(output1j$Px1)) {
  tmp1<-myfunctionE3(0.2,0.2,output1j$Px1[i],0.2,0.5,1,output1j$k[i],2,2,0.1,10000,10000,output1j$siglevel[i],1)
  pow<-tmp1$power2
  d1[i]<-pow
}
power2<-as.data.frame(do.call(rbind, d1))
output1j$power2 <- power2$V1

#draw fig6A
data <- output1j[1:63,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,21,1.5),seq(23,41,1.5),seq(44,62,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("k", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("k", "Px1"),
                                     names_to = "group", values_to = "power")  

p1<-ggplot()+
  geom_line(data=data_long,aes(x = k, y = power, linetype = group, color = Px1),linewidth=1)+
  geom_point(data=data_long2,aes(x = k, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("k")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,2),breaks = seq(0,2,0.2))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#draw fig6B
data <- output1j[64:126,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,21,1.5),seq(23,41,1.5),seq(44,62,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("k", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("k", "Px1"),
                                     names_to = "group", values_to = "power")  

p2<-ggplot()+
  geom_line(data=data_long,aes(x = k, y = power, linetype = group, color = Px1),linewidth=1)+
  geom_point(data=data_long2,aes(x = k, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("k")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,2),breaks = seq(0,2,0.2))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#draw fig6C
data <- output1j[127:189,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,21,1.5),seq(23,41,1.5),seq(44,62,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("k", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("k", "Px1"),
                                     names_to = "group", values_to = "power")  

p3<-ggplot()+
  geom_line(data=data_long,aes(x = k, y = power, linetype = group, color = Px1),linewidth=1)+
  geom_point(data=data_long2,aes(x = k, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("k")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,2),breaks = seq(0,2,0.2))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#draw fig6D
data <- output1j[190:252,-3]
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,21,1.5),seq(23,41,1.5),seq(44,62,1.5)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("k", "Px1"),
                                   names_to = "group", values_to = "power")  
data_long2 <- data2 %>% pivot_longer(cols = -c("k", "Px1"),
                                     names_to = "group", values_to = "power")  

p4<-ggplot()+
  geom_line(data=data_long,aes(x = k, y = power, linetype = group, color = Px1),linewidth=1)+
  geom_point(data=data_long2,aes(x = k, y = power, color=Px1, shape=Px1),size=3)+
  scale_shape_manual(values = c(15, 19, 17))+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("k")+
  ylab("power")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,2),breaks = seq(0,2,0.2))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))

#save fig.6A, fig.6B, fig.6C, and fig.6D
ggsave("power,power2&k-px1-1.pdf", plot = p1, height = 4, width = 5)
ggsave("power,power2&k-px1-2.pdf", plot = p2, height = 4, width = 5)
ggsave("power,power2&k-px1-3.pdf", plot = p3, height = 4, width = 5)
ggsave("power,power2&k-px1-4.pdf", plot = p4, height = 4, width = 5)





