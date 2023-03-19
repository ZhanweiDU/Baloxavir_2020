withinHostTrans = function(tGap, tGap6days, epsilonSet, betaScaling, pScaling, VsScaling,  
                           cScaling, ZScaling, deltaScaling, rScaling, sigmaScaling, HourScaling, modelFlag){
  
  HoursDay = 24*2
  if(!missing(modelFlag)){
    if( length(grep(modelFlag, 'withinHost'))>0 ){
      tPeriod = (HoursDay*(10)*HourScaling)
    }else{
      tPeriod= (HoursDay*(20)*HourScaling)
    }
  }else{
    tPeriod= (HoursDay*(20)*HourScaling)
  }
  
  Is = 0; Ir = 0; Vr = 0; T = 4*10^8;
  
  mu = 0
  sigma= 0.1*sigmaScaling;
  epsilon = 0;
  
  beta = ( 9.9*10^(-2))/HoursDay/HourScaling*10^(betaScaling) ; # b Infection rate, can adjust the peak timing
  p = (1.2*10^(-5))/HoursDay/HourScaling *10^(pScaling); #p Virus production rate, can adjust the peak value
  Vs= (7.7*10^(-3)*10^(VsScaling));  # Initial sensitive viral load
  c =  ( 8.1*10^(-2))/HoursDay/HourScaling*10^(cScaling); # dV
  Z= 0* (3.4*10^(-1))*10^(ZScaling);# X0
  delta =  ( 5.0*10^(-1) )/HoursDay/HourScaling*10^(deltaScaling);#dI
  r = 0* 1/HoursDay/HourScaling*10^(rScaling);#10^(rScaling); # r Immune response growth rate
  
  epsilon = 0;##a
  k = 1/HoursDay/HourScaling
  
  for (t in 1:tPeriod){## 20 minutes over 7 days
    if(t==tGap*HourScaling+1 & epsilonSet>0){
      epsilon = epsilonSet;
    }
  
    dT  = -beta*T[t]*(Vs[t]+Vr[t]);
    dIs = beta*T[t]*Vs[t] - delta*Is[t];
    dIr = beta*T[t]*Vr[t] - delta*Ir[t];
    dZ  = r*Z[t];
    dVs = (1-epsilon)*(1-mu)*p*Is[t] - c*Vs[t] - k*Z[t]*Vs[t];
    dVr = (1-epsilon)*mu*p*Is[t] + (1-sigma)*p*Ir[t] - c*Vr[t] - k*Z[t]*Vr[t];
    
    T[t+1]  = T[t] + dT;
    Is[t+1] = Is[t]+ dIs;
    Ir[t+1] = Ir[t]+ dIr;
    Z[t+1]  = Z[t] + dZ;
    Vs[t+1] = Vs[t]+ dVs;
    Vr[t+1] = Vr[t]+ dVr;
    
    
    if (T[t+1]<0)  {T[t+1]=0; }
    if (Is[t+1]<0) {Is[t+1]=0; }
    if (Ir[t+1]<0) {Ir[t+1]=0; }
    if (Z[t+1]<0)  {Z[t+1]=0; }
    if (Vs[t+1]>Vs[t] & t>48) { Vs[t+1]=Vs[t]; }
    if (Vs[t+1]<0) {Vs[t+1]=0; }
    if (Vr[t+1]<0) {Vr[t+1]=0; }
    
  }
  Vrs = Vr+Vs;
  Vrs = Vrs[ seq(1, tGap6days*HoursDay/24, by=HoursDay/24)];
}
