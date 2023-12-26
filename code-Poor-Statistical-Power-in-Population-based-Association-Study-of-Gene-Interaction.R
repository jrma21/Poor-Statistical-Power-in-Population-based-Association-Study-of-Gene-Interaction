#basic function
myfunctionE<-function(Pm1,Pm2,Px1,Px2,d1,d2,k,r1,r2,PD,nD,nd,sig.level,df){   
  #Calculate P(x1=1│D) and P(x2=1│D)
  Px1d<-r1*Px1/(1+(r1-1)*Px1)           
  Px2d<-r2*Px2/(1+(r2-1)*Px2)           
  #Calculate P(x1=1, x2=1│D),P(x1=1, x2=0│D),P(x1=0, x2=1│D),and P(x1=0, x2=0│D)
  P11d<-(1+k)*Px1d*Px2d                
  P10d<-Px1d-(1+k)*Px1d*Px2d            
  P01d<-Px2d-(1+k)*Px1d*Px2d            
  P00d<-1-Px1d-Px2d+(1+k)*Px1d*Px2d  
  #Calculate t1,t2 (t1,t2 are parameters required to calculate conditional probabilities in the control group)
  ta<-1-Px1d*PD/Px1                    
  tb<-1-(1-Px1d)*PD/(1-Px1)
  tc<-1-Px2d*PD/Px2
  td<-1-(1-Px2d)*PD/(1-Px2)
  t1<-ta/tb
  t2<-tc/td
  #Calculate P(x1=1│d) and P(x2=1│d)
  dPx1<-t1*Px1/(1+(t1-1)*Px1)          
  dPx2<-t2*Px2/(1+(t2-1)*Px2)    
  #Calculate P(x1=1,x2=1│d),P(x1=1,x2=0│d),P(x1=0,x2=1│d),and P(x1=0, x2=0│d)
  dP11<-(1+k)*dPx1*dPx2                 
  dP10<-dPx1-(1+k)*dPx1*dPx2            
  dP01<-dPx2-(1+k)*dPx1*dPx2           
  dP00<-1-dPx1-dPx2+(1+k)*dPx1*dPx2  
  #Calculate P(m1=1|x1=1),P(m2=1|x2=1),P(m1=1|x1=0),P(m2=1|x2=0),P(m1=0|x1=1),P(m2=0|x2=1),P(m1=0|x1=0), and P(m2=0|x2=0)
  Pa1<-(Pm1*Px1+d1)/Px1                
  Pa2<-(Pm2*Px2+d2)/Px2                 
  Pb1<-(Pm1*(1-Px1)+d1)/(1-Px1)         
  Pb2<-(Pm2*(1-Px2)+d2)/(1-Px2)         
  Pc1<-((1-Pm1)*Px1+d1)/Px1             
  Pc2<-((1-Pm2)*Px2+d2)/Px2            
  Pd1<-((1-Pm1)*(1-Px1)+d1)/(1-Px1)    
  Pd2<-((1-Pm2)*(1-Px2)+d2)/(1-Px2)    
  #Calculate 4 frequencies f11D,f22D,f12D, and f21D for the case group
  f11D<-Pa1*Pa2*P11d+Pb1*Pa2*P01d+Pa1*Pb2*P10d+Pb1*Pb2*P00d 
  f22D<-Pc1*Pc2*P11d+Pd1*Pc2*P01d+Pc1*Pd2*P10d+Pd1*Pd2*P00d
  f12D<-Pa1*Pc2*P11d+Pb1*Pc2*P01d+Pa1*Pd2*P10d+Pb1*Pd2*P00d
  f21D<-Pc1*Pa2*P11d+Pd1*Pa2*P01d+Pc1*Pb2*P10d+Pd1*Pb2*P00d
  #Calculate 4 frequencies f11d,f22d,f12d, and f21d for the control group
  f11d<-Pa1*Pa2*dP11+Pb1*Pa2*dP01+Pa1*Pb2*dP10+Pb1*Pb2*dP00 
  f22d<-Pc1*Pc2*dP11+Pd1*Pc2*dP01+Pc1*Pd2*dP10+Pd1*Pd2*dP00
  f12d<-Pa1*Pc2*dP11+Pb1*Pc2*dP01+Pa1*Pd2*dP10+Pb1*Pd2*dP00
  f21d<-Pc1*Pa2*dP11+Pd1*Pa2*dP01+Pc1*Pb2*dP10+Pd1*Pb2*dP00
  #Calculate the expectation E1 of the difference of frequencies for the case group and the expectation E2 of the difference of frequencies for the control group 
  E1<-f11D*f22D-f12D*f21D              
  E2<-f11d*f22d-f12d*f21d    
  #Calculate the expectation dE of the effective size
  dE<-E1-E2 
  #Calculate the covariance V1 of the difference of frequencies for the case group and the covariance V2 of the difference of frequencies for the control group
  V1<-((f11D+f12D)*(1-f11D-f12D)*(f11D+f21D)*(1-f11D-f21D)+(1-2*f11D-2*f12D)*(1-2*f11D-2*f21D)*(f11D*f22D-f12D*f21D)-(f11D*f22D-f12D*f21D)*(f11D*f22D-f12D*f21D))/nD 
  V2<-((f11d+f12d)*(1-f11d-f12d)*(f11d+f21d)*(1-f11d-f21d)+(1-2*f11d-2*f12d)*(1-2*f11d-2*f21d)*(f11d*f22d-f12d*f21d)-(f11d*f22d-f12d*f21d)*(f11d*f22d-f12d*f21d))/nd     
  #Calculate V, the sum of V1 and V2
  V<-V1+V2  
  #Calculate the statistic T
  T<-dE*dE/V
  #Convert the statistic into effect size W (k1 is the minimum value in the number of rows and columns of the contingency table, k1 is 2 and K1-1 is 1, so they are omitted from the denominator)
  W<-sqrt(T/(nD+nd))    
  #Connect sig.level, df, and power1
  k1<- qchisq(sig.level, df = df, lower = FALSE)
  #Calculate the power1 value in the chi-square test
  power<-pchisq(k1, df = df, ncp = (nd+nD)*W^2, lower = FALSE)
  #Output E1,E2,dE,V1,V2,V,T,W, and power1
  cat("E1=",E1,"\n","E2=",E2,"\n","dE=",dE,"\n","V1=",V1,"\n","V2=",V2,"\n","V=",V,"\n","T=",T,"\n","W=",W,"\n","power1=",power1) 
}

myfunctionE2<-function(Pm1,Pm2,Px1,Px2,d1a,d2a,k,r1,r2,PD,nD,nd,sig.level,df){   #d1a=D1*D1,d2a=D2*D2
  #Calculate d1max and d2max
  d1max<-min(Px1*(1-Pm1),Pm1*(1-Px1))           
  d2max<-min(Px2*(1-Pm2),Pm2*(1-Px2))
  #Calculate d1 and d2
  d1<-d1max*d1a          
  d2<-d2max*d2a           
  #Calculate r1a and r2a (r1a=r1'*r1',r2a=r2'*r2')
  r1a<-d1*d1/((Px1*(1-Pm1))*Pm1*(1-Px1))                
  r2a<-d2*d2/((Px2*(1-Pm2))*Pm2*(1-Px2)) 
  #Calculate P(x1=1│D) and P(x2=1│D)
  Px1d<-r1*Px1/(1+(r1-1)*Px1)           
  Px2d<-r2*Px2/(1+(r2-1)*Px2)           
  #Calculate P(x1=1, x2=1│D),P(x1=1, x2=0│D),P(x1=0, x2=1│D),and P(x1=0, x2=0│D)
  P11d<-(1+k)*Px1d*Px2d                
  P10d<-Px1d-(1+k)*Px1d*Px2d            
  P01d<-Px2d-(1+k)*Px1d*Px2d            
  P00d<-1-Px1d-Px2d+(1+k)*Px1d*Px2d  
  #Calculate t1,t2 (t1,t2 are parameters required to calculate conditional probabilities in the control group)
  ta<-1-Px1d*PD/Px1                    
  tb<-1-(1-Px1d)*PD/(1-Px1)
  tc<-1-Px2d*PD/Px2
  td<-1-(1-Px2d)*PD/(1-Px2)
  t1<-ta/tb
  t2<-tc/td
  #Calculate P(x1=1│d) and P(x2=1│d)
  dPx1<-t1*Px1/(1+(t1-1)*Px1)          
  dPx2<-t2*Px2/(1+(t2-1)*Px2)    
  #Calculate P(x1=1,x2=1│d),P(x1=1,x2=0│d),P(x1=0,x2=1│d),and P(x1=0, x2=0│d)
  dP11<-(1+k)*dPx1*dPx2                 
  dP10<-dPx1-(1+k)*dPx1*dPx2            
  dP01<-dPx2-(1+k)*dPx1*dPx2           
  dP00<-1-dPx1-dPx2+(1+k)*dPx1*dPx2  
  #Calculate P(m1=1|x1=1),P(m2=1|x2=1),P(m1=1|x1=0),P(m2=1|x2=0),P(m1=0|x1=1),P(m2=0|x2=1),P(m1=0|x1=0), and P(m2=0|x2=0)
  Pa1<-(Pm1*Px1+d1)/Px1                
  Pa2<-(Pm2*Px2+d2)/Px2                 
  Pb1<-(Pm1*(1-Px1)+d1)/(1-Px1)         
  Pb2<-(Pm2*(1-Px2)+d2)/(1-Px2)         
  Pc1<-((1-Pm1)*Px1+d1)/Px1             
  Pc2<-((1-Pm2)*Px2+d2)/Px2            
  Pd1<-((1-Pm1)*(1-Px1)+d1)/(1-Px1)    
  Pd2<-((1-Pm2)*(1-Px2)+d2)/(1-Px2)    
  #Calculate 4 frequencies f11D,f22D,f12D, and f21D for the case group
  f11D<-Pa1*Pa2*P11d+Pb1*Pa2*P01d+Pa1*Pb2*P10d+Pb1*Pb2*P00d 
  f22D<-Pc1*Pc2*P11d+Pd1*Pc2*P01d+Pc1*Pd2*P10d+Pd1*Pd2*P00d
  f12D<-Pa1*Pc2*P11d+Pb1*Pc2*P01d+Pa1*Pd2*P10d+Pb1*Pd2*P00d
  f21D<-Pc1*Pa2*P11d+Pd1*Pa2*P01d+Pc1*Pb2*P10d+Pd1*Pb2*P00d
  #Calculate 4 frequencies f11d,f22d,f12d, and f21d for the control group
  f11d<-Pa1*Pa2*dP11+Pb1*Pa2*dP01+Pa1*Pb2*dP10+Pb1*Pb2*dP00 
  f22d<-Pc1*Pc2*dP11+Pd1*Pc2*dP01+Pc1*Pd2*dP10+Pd1*Pd2*dP00
  f12d<-Pa1*Pc2*dP11+Pb1*Pc2*dP01+Pa1*Pd2*dP10+Pb1*Pd2*dP00
  f21d<-Pc1*Pa2*dP11+Pd1*Pa2*dP01+Pc1*Pb2*dP10+Pd1*Pb2*dP00
  #Calculate the expectation E1 of the difference of frequencies for the case group and the expectation E2 of the difference of frequencies for the control group 
  E1<-f11D*f22D-f12D*f21D              
  E2<-f11d*f22d-f12d*f21d    
  #Calculate the expectation dE of the effective size
  dE<-E1-E2 
  #Calculate the covariance V1 of the difference of frequencies for the case group and the covariance V2 of the difference of frequencies for the control group
  V1<-((f11D+f12D)*(1-f11D-f12D)*(f11D+f21D)*(1-f11D-f21D)+(1-2*f11D-2*f12D)*(1-2*f11D-2*f21D)*(f11D*f22D-f12D*f21D)-(f11D*f22D-f12D*f21D)*(f11D*f22D-f12D*f21D))/nD 
  V2<-((f11d+f12d)*(1-f11d-f12d)*(f11d+f21d)*(1-f11d-f21d)+(1-2*f11d-2*f12d)*(1-2*f11d-2*f21d)*(f11d*f22d-f12d*f21d)-(f11d*f22d-f12d*f21d)*(f11d*f22d-f12d*f21d))/nd     
  #Calculate V, the sum of V1 and V2
  V<-V1+V2  
  #Calculate the statistic T
  T<-dE*dE/V
  #Convert the statistic into effect size W (k1 is the minimum value in the number of rows and columns of the contingency table, k1 is 2 and K1-1 is 1, so they are omitted from the denominator)
  W<-sqrt(T/(nD+nd))    
  #Connect sig.level, df, and powe1r
  k1<- qchisq(sig.level, df = df, lower = FALSE)
  #Calculate the power value in the chi-square test
  power<-pchisq(k1, df = df, ncp = (nd+nD)*W^2, lower = FALSE)
  #Output d1,d2,r1a,r2a,and power
  cat("d1=",d1,"\n","d2=",d2,"\n","r1a=",r1a,"\n","r2a=",r2a,"power1=",power1) 
}

myfunctionE3<-function(Pm1,Pm2,Px1,Px2,d1a,d2a,k,r1,r2,PD,nD,nd,sig.level,df){   #d1a=D1*D1,d2a=D2*D2
  #Calculate d1max and d2max
  d1max<-min(Px1*(1-Pm1),Pm1*(1-Px1))           
  d2max<-min(Px2*(1-Pm2),Pm2*(1-Px2))
  #Calculate d1 and d2
  d1<-d1max*d1a          
  d2<-d2max*d2a           
  #Calculate r1a and r2a (r1a=r1'*r1',r2a=r2'*r2')
  r1a<-d1*d1/((Px1*(1-Pm1))*Pm1*(1-Px1))                
  r2a<-d2*d2/((Px2*(1-Pm2))*Pm2*(1-Px2)) 
  #Calculate P(x1=1│D) and P(x2=1│D)
  Px1d<-r1*Px1/(1+(r1-1)*Px1)           
  Px2d<-r2*Px2/(1+(r2-1)*Px2)           
  #Calculate P(x1=1, x2=1│D),P(x1=1, x2=0│D),P(x1=0, x2=1│D),and P(x1=0, x2=0│D)
  P11d<-(1+k)*Px1d*Px2d                
  P10d<-Px1d-(1+k)*Px1d*Px2d            
  P01d<-Px2d-(1+k)*Px1d*Px2d            
  P00d<-1-Px1d-Px2d+(1+k)*Px1d*Px2d  
  #Calculate t1,t2 (t1,t2 are parameters required to calculate conditional probabilities in the control group)
  ta<-1-Px1d*PD/Px1                    
  tb<-1-(1-Px1d)*PD/(1-Px1)
  tc<-1-Px2d*PD/Px2
  td<-1-(1-Px2d)*PD/(1-Px2)
  t1<-ta/tb
  t2<-tc/td
  #Calculate P(x1=1│d) and P(x2=1│d)
  dPx1<-t1*Px1/(1+(t1-1)*Px1)          
  dPx2<-t2*Px2/(1+(t2-1)*Px2)    
  #Calculate P(x1=1,x2=1│d),P(x1=1,x2=0│d),P(x1=0,x2=1│d),and P(x1=0, x2=0│d)
  dP11<-(1+k)*dPx1*dPx2                 
  dP10<-dPx1-(1+k)*dPx1*dPx2            
  dP01<-dPx2-(1+k)*dPx1*dPx2           
  dP00<-1-dPx1-dPx2+(1+k)*dPx1*dPx2  
  #Calculate P(m1=1|x1=1),P(m2=1|x2=1),P(m1=1|x1=0),P(m2=1|x2=0),P(m1=0|x1=1),P(m2=0|x2=1),P(m1=0|x1=0), and P(m2=0|x2=0)
  Pa1<-(Pm1*Px1+d1)/Px1                
  Pa2<-(Pm2*Px2+d2)/Px2                 
  Pb1<-(Pm1*(1-Px1)+d1)/(1-Px1)         
  Pb2<-(Pm2*(1-Px2)+d2)/(1-Px2)         
  Pc1<-((1-Pm1)*Px1+d1)/Px1             
  Pc2<-((1-Pm2)*Px2+d2)/Px2            
  Pd1<-((1-Pm1)*(1-Px1)+d1)/(1-Px1)    
  Pd2<-((1-Pm2)*(1-Px2)+d2)/(1-Px2)    
  #Calculate 4 frequencies f11D,f22D,f12D, and f21D for the case group
  f11D<-Pa1*Pa2*P11d+Pb1*Pa2*P01d+Pa1*Pb2*P10d+Pb1*Pb2*P00d 
  f22D<-Pc1*Pc2*P11d+Pd1*Pc2*P01d+Pc1*Pd2*P10d+Pd1*Pd2*P00d
  f12D<-Pa1*Pc2*P11d+Pb1*Pc2*P01d+Pa1*Pd2*P10d+Pb1*Pd2*P00d
  f21D<-Pc1*Pa2*P11d+Pd1*Pa2*P01d+Pc1*Pb2*P10d+Pd1*Pb2*P00d
  #Calculate 4 frequencies f11d,f22d,f12d, and f21d for the control group
  f11d<-Pa1*Pa2*dP11+Pb1*Pa2*dP01+Pa1*Pb2*dP10+Pb1*Pb2*dP00 
  f22d<-Pc1*Pc2*dP11+Pd1*Pc2*dP01+Pc1*Pd2*dP10+Pd1*Pd2*dP00
  f12d<-Pa1*Pc2*dP11+Pb1*Pc2*dP01+Pa1*Pd2*dP10+Pb1*Pd2*dP00
  f21d<-Pc1*Pa2*dP11+Pd1*Pa2*dP01+Pc1*Pb2*dP10+Pd1*Pb2*dP00
  #Calculate the expectation E3 of the odd ratio for the case group and the control group 
  E3<-(f11D+f12D)*(f21d+f22d)/((f21D+f22D)*(f11d+f12d))    
  #Calculate the variance of the estimated variance of the logarithm of the odd ratio
  V3<-1/(nD*(f11D+f12D))+1/(nD*(f21D+f22D))+1/(nd*(f11d+f12d))+1/(nd*(f21d+f22d))
  #Calculate the statistic T2
  T2<-log(E3)*log(E3)/V3
  #Convert the statistic into effect size W (k1 is the minimum value in the number of rows and columns of the contingency table, k1 is 2 and K1-1 is 1, so they are omitted from the denominator)
  W2<-sqrt(T2/(nD+nd))    
  #Connect sig.level, df, and power2
  k1<- qchisq(sig.level, df = df, lower = FALSE)
  #Calculate the power value in the chi-square test
  power2<-pchisq(k1, df = df, ncp = (nd+nD)*W2^2, lower = FALSE)
  #Output d1,d2,r1a,r2a,and power2
  cat("d1=",d1,"\n","d2=",d2,"\n","r1a=",r1a,"\n","r2a=",r2a,"power2=",power2,"k=",k) 
}

#data of fig.1
#investigate the relationship between variance（V）and the total chromosome number of the disease population (nD) under different values of the frequency of allele X1 (Px1)  

z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,2,2,0.1,seq(3450,50000,2450),0.05,1) 
z <- myfunctionE(0.05,0.05,0.1,0.05,0.0475,0.0475,1,2,2,0.1,seq(3450,50000,2450),0.05,1)
z <- myfunctionE(0.05,0.05,0.15,0.05,0.0475,0.0475,1,2,2,0.1,seq(3450,50000,2450),0.05,1)

#data of fig.2
#investigate the relationship between mathematical expectation (dE) , the relative risk ratio of X1 (r1,r2) under different values of the gene interaction strength of X1 and X2 (k)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,2,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,0.05,1) 
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,0.05,1) 
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,0.5,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,0.05,1) 

#investigate the relationship between dE and k under different values of r1,r2(r1=r2)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2.5,2.5,0.1,10000,0.05,1) 
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2,2,0.1,10000,0.05,1) 
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),1.5,1.5,0.1,10000,0.05,1) 

#data of fig.3

#investigate the relationship between power1 and r1,r2(r1=r2) under different values of k(fig.3A)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,2,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,10000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,10000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,0.5,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,10000,0.05,1)

#investigate the relationship between power1 and r1,r2(r1=r2) under different values of nD(fig.3B)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,seq(1,4,0.15),seq(1,4,0.15),0.1,10000,10000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,seq(1,4,0.15),seq(1,4,0.15),0.1,5000,5000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,1,seq(1,4,0.15),seq(1,4,0.15),0.1,1000,1000,0.05,1)

#data of fig.4

#investigate the relationship between power1 and k under different values of r1,r2(r1=r2)(fig.4A)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2.5,2.5,0.1,10000,10000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2,2,0.1,10000,10000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),1.2,1.2,0.1,10000,10000,0.05,1)

#investigate the relationship between power1 and k under different values of nD(fig.4B)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2,2,0.1,10000,10000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2,2,0.1,5000,5000,0.05,1)
z <- myfunctionE(0.05,0.05,0.05,0.05,0.0475,0.0475,seq(0,5,0.25),2,2,0.1,1000,1000,0.05,1)

#data of fig.5

#investigate the relationship between the linkage disequilibrium value of M1 and X1 (d1), r1a, and under different values of Px1

myfunctionE2(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)

#data of fig.6

#investigate the relationship between power1,power2 and d1a (d1/d1max) under different values of Px1 and level of significance (sig.level)

#fig6A

myfunctionE3(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE3(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE3(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.05,1)

#fig6B

myfunctionE3(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.0005,1)
myfunctionE3(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.0005,1)
myfunctionE3(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.0005,1)
myfunctionE2(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.0005,1)
myfunctionE2(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.0005,1)
myfunctionE2(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.0005,1)

#fig6C

myfunctionE3(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.000005,1)
myfunctionE3(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.000005,1)
myfunctionE3(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.000005,1)
myfunctionE2(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.000005,1)
myfunctionE2(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.000005,1)
myfunctionE2(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.000005,1)

#fig6D

myfunctionE3(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.00000005,1)
myfunctionE3(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.00000005,1)
myfunctionE3(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.00000005,1)
myfunctionE2(0.2,0.2,0.2,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.00000005,1)
myfunctionE2(0.2,0.2,0.3,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.00000005,1)
myfunctionE2(0.2,0.2,0.1,0.2,seq(0.05,1,0.05),1,1,2,2,0.1,10000,10000,0.00000005,1)

#data of fig.7

#investigate the relationship between power1,power2 and d1a under different values of Px1 and k

#fig7A

myfunctionE3(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.05,1)
myfunctionE3(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.05,1)
myfunctionE3(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.05,1)
myfunctionE2(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.05,1)

#fig7B

myfunctionE3(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.0005,1)
myfunctionE3(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.0005,1)
myfunctionE3(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.0005,1)
myfunctionE2(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.0005,1)
myfunctionE2(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.0005,1)
myfunctionE2(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.0005,1)

#fig7C

myfunctionE3(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.000005,1)
myfunctionE3(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.000005,1)
myfunctionE3(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.000005,1)
myfunctionE2(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.000005,1)
myfunctionE2(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.000005,1)
myfunctionE2(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.000005,1)

#fig7D

myfunctionE3(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.00000005,1)
myfunctionE3(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.00000005,1)
myfunctionE3(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.00000005,1)
myfunctionE2(0.2,0.2,0.2,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.00000005,1)
myfunctionE2(0.2,0.2,0.3,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.00000005,1)
myfunctionE2(0.2,0.2,0.1,0.2,0.5,0.5,seq(0,2,0.1),2,2,0.1,10000,10000,0.00000005,1)

setwd("")
library(reshape2) 
library(ggplot2)   
library(ggpubr)

#fig1.investigate the relationship between variance（V）and the total chromosome number of the disease population (nD) under different values of the frequency of allele X1 (Px1) 

data<-read.table("clipboard",head=T)
data$Px1<-as.factor(data$Px1)
colnames(data)[4]<-c("mlogv")
expv<-exp(data$mlogv)
data[,5]<-expv
colnames(data)[5]<-c("expmlogv")
data2<-data[-seq(2,59,2),]

p1<-ggplot()+
  geom_line(data = data,aes(x=n,y=v,group=Px1,color=Px1,shape=Px1),linewidth=0.5)+
  geom_point(data = data2,aes(x=n,y=v,group=Px1,color=Px1,shape=Px1),size=2)+
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
  scale_x_continuous(expand=c(0,0),limits = c(0,52000),breaks = seq(0,50000,5000))

ggsave("power-n-px1.pdf", plot = p1, height = 3, width = 4)
write.csv(data,"variance-n.csv")

#fig2

#investigate the relationship between mathematical expectation (dE) , the relative risk ratio of X1 (r1,r2) under different values of the gene interaction strength of X1 and X2 (k) (fig.2A)

data<-read.table("clipboard",head=T)
data$RR<-as.factor(data$RR)
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=K,y=size,group=RR,color=RR,shape=RR),linewidth=1.3)+
  geom_point(data = data2,aes(x=K,y=size,group=RR,color=RR,shape=RR),size=3)+
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

#investigate the relationship between dE and k under different values of r1,r2(r1=r2) (fig.2B)
data<-read.table("clipboard",head=T)
data$K<-as.factor(data$K)
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=R,y=size,group=K,color=K,shape=K),linewidth=1.3)+
  geom_point(data = data2,aes(x=R,y=size,group=K,color=K,shape=K),size=3)+
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("RR")+
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



#fig3.

#investigate the relationship between power1 and r1,r2(r1=r2) under different values of k(fig.3A)
data<-read.table("clipboard",head=T)
data$k<-factor(data$k,levels=c("2","1","0.5"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=R,y=power,group=k,color=k,shape=k),linewidth=1.3)+
  geom_point(data = data2,aes(x=R,y=power,group=k,color=k,shape=k),size=3)+ 
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
data<-read.table("clipboard",head=T)
data$n<-factor(data$n,levels=c("10000","5000","1000"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=R,y=power,group=n,color=n,shape=n),linewidth=1.3)+
  geom_point(data = data2,aes(x=R,y=power,group=n,color=n,shape=n),size=3)+ 
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



#fig4.power-r

#investigate the relationship between power1 and k under different values of r1,r2(r1=r2)(fig.4A)
data<-read.table("clipboard",head=T)
data$r<-factor(data$r,levels=c("2.5","2","1.2"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=K,y=power,group=r,color=r,shape=r),linewidth=1.3)+
  geom_point(data = data2,aes(x=K,y=power,group=r,color=r,shape=r),size=3)+ 
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
data<-read.table("clipboard",head=T)
data$n<-factor(data$n,levels=c("10000","5000","1000"))
data2<-data[-c(1,22,43,seq(2,20,2),seq(23,41,2),seq(44,62,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=K,y=power,group=n,color=n,shape=n),linewidth=1.3)+
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



#fig5.power-r

#investigate the relationship between the linkage disequilibrium value of M1 and X1 (d1), r1a, and under different values of Px1
#d1a-d1-Px1(fig5A)
data<-read.table("clipboard",head=T)
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,2),seq(22,38,2),seq(42,58,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=d1a,y=d1,group=Px1,color=Px1,shape=Px1),linewidth=1.3)+
  geom_point(data = data2,aes(x=d1a,y=d1,group=Px1,color=Px1,shape=Px1),size=3)+ 
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("d1a")+
  ylab("d1")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0,1.05),breaks = seq(0,1.05,0.1))+
  scale_y_continuous(limits = c(0,0.17),breaks = seq(0,0.17,0.05))
ggsave("d1a-d1-Px1.pdf", plot = p1, height = 4, width = 5)

#d1a-r1a-px1
data<-read.table("clipboard",head=T)
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,2),seq(22,38,2),seq(42,58,2)),]

p1<-ggplot()+
  geom_line(data = data,aes(x=d1a,y=r1a,group=Px1,color=Px1,shape=Px1),linewidth=1.3)+
  geom_point(data = data2,aes(x=d1a,y=r1a,group=Px1,color=Px1,shape=Px1),size=3)+ 
  scale_color_manual(values = c('#103666','#00ae71','#f5b06d'))+ 
  xlab("d1a")+
  ylab("r1a")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.box.background = element_rect(color="black"),
        legend.key.size=unit(20,"pt"))+
  scale_x_continuous(expand=c(0,0),limits = c(0,1.05),breaks = seq(0,1.05,0.1))+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))
ggsave("d1a-r1a-px1.pdf", plot = p1, height = 4, width = 5)

#fig6.investigate the relationship between power1,power2 and d1a (d1/d1max) under different values of Px1 and level of significance (sig.level)
data<-read.table("clipboard",head=T)
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,2),seq(22,38,2),seq(42,58,2)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("d1a", "Px1"),
                                   names_to = "group", values_to = "power")  

p1<-ggplot(data_long)+
  geom_line(aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = d1a, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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

p2<-ggplot(data_long)+
  geom_line(aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = d1a, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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

p3<-ggplot(data_long)+
  geom_line(aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = d1a, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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

p4<-ggplot(data_long)+
  geom_line(aes(x = d1a, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = d1a, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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

ggsave("power,power2&dla-px1-1.pdf", plot = p1, height = 4, width = 5)
ggsave("power,power2&dla-px1-2.pdf", plot = p2, height = 4, width = 5)
ggsave("power,power2&dla-px1-3.pdf", plot = p3, height = 4, width = 5)
ggsave("power,power2&dla-px1-4.pdf", plot = p4, height = 4, width = 5)

#fig7.investigate the relationship between power1,power2 and d1a under different values of Px1 and k
data<-read.table("clipboard",head=T)
data$Px1<-factor(data$Px1,levels=c("0.2","0.3","0.1"))
data2<-data[-c(seq(2,18,2),seq(22,38,2),seq(42,58,2)),]

library(tidyverse)
data_long <- data %>% pivot_longer(cols = -c("k", "Px1"),
                                   names_to = "group", values_to = "power")  

p1<-ggplot(data_long)+
  geom_line(aes(x = k, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = k, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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

p2<-ggplot(data_long)+
  geom_line(aes(x = k, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = k, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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

p3<-ggplot(data_long)+
  geom_line(aes(x = k, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = k, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.1))

p4<-ggplot(data_long)+
  geom_line(aes(x = k, y = power, linetype = group, color = Px1),linewidth=0.8)+
  geom_point(aes(x = k, y = power, linetype = group, color=Px1, shape=Px1),size=2)+
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
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.1))

ggsave("power,power2&k-px1-1.pdf", plot = p1, height = 4, width = 5)
ggsave("power,power2&k-px1-2.pdf", plot = p2, height = 4, width = 5)
ggsave("power,power2&k-px1-3.pdf", plot = p3, height = 4, width = 5)
ggsave("power,power2&k-px1-4.pdf", plot = p4, height = 4, width = 5)



















