## Description:
# Agent-based simulations
# (1) Simualte each agent, with (1) random time for symptom since infection and (2) random time for treatment since symptom
# (2) Estimate MCMC likelihood as Joint Normal probability of agents with observed mean and stds.

optimal_withinHostTrans_triple_Median_Symptom = function(x, PopDist_Baloxavir, BaloxavirScaleHour,HourScaling,
                                                         PopDist_Placebo, PlaceboScaleHour, PopDist_Oseltamivir, 
                                                         OseltamivirScaleHour, lb, ub, returnLog, eiSingle, parasNum, otherParas, outputTable){
  # set.seed( sum(x))
  # print(x)
  sigmaScaling = 1;
  sdScaling = 1;
  if(parasNum>=7){
    betaScaling = x[1];
    pScaling = x[2]; # 3 times
    VsScaling = x[3]; # 3 times
    cScaling = x[4];
    ZScaling = x[5];
    deltaScaling = x[6];
    rScaling = x[7];
    SymptomScaling = c(0, 0, 0);
    epsilonSets = c(0)/10;
    
    FlagLim = FALSE
  }else{
    betaScaling = otherParas[1];
    pScaling = otherParas[2]; # 3 times
    VsScaling = otherParas[3]; # 3 times
    cScaling = otherParas[4];
    ZScaling = otherParas[5];
    deltaScaling = otherParas[6];
    rScaling = otherParas[7];
    SymptomScaling = c(0, 0, 0);
    epsilonSets = c(0, x[1], x[1])/10;
    
    FlagLim = x[1]< lb[1] || x[1]>ub[1]
  }
  
  if(FlagLim)
  {
    return(-10^10)
  }
  else{
    
    
    logliks = 0;
    
    for (ei in eiSingle:eiSingle){
      if(ei==1){
        temp_PlaceboScaleHour = c(PlaceboScaleHour, rep(tail(PlaceboScaleHour,1), 24*3 ) );
        PopDist = PopDist_Placebo;
      }   
      if(ei==2){
        temp_PlaceboScaleHour = c(BaloxavirScaleHour, rep(tail(BaloxavirScaleHour,1), 24*3 ) );
        PopDist = PopDist_Baloxavir;
      }
      if(ei==3){
        temp_PlaceboScaleHour = c(OseltamivirScaleHour, rep(tail(OseltamivirScaleHour,1), 24*3 ) );
        PopDist = PopDist_Oseltamivir;
      }
      
      epsilonSet = epsilonSets[ei];
      temp_SymptomScaling = SymptomScaling[ei];
      Vrs_all = c();
      tGapAll = c()
      Vrs_allPeriod = c()
      for(i in 1:length(PopDist)){
        tempPop = PopDist[i];
        for(j in 1:tempPop){
          tempA =  round(rgamma(1, shape=3.9855, rate = 1/6.3190))
          tGap = tempA + round(rlnorm(1, meanlog=0.3686, sdlog = 0.4052*sdScaling)*24);# + round(48*temp_SymptomScaling);
          tGapAll = c(tGapAll, tempA)
          
          tGap6days = tGap+9*24; #% 5 additional days
          
          HourScaling=10
          Vrs = withinHostTransRCpp240(tGap, tGap6days, epsilonSet, betaScaling, pScaling, VsScaling,
                                       cScaling, ZScaling, deltaScaling, rScaling, sigmaScaling, HourScaling); #% tGap=0 means the starting usage of drug
          
          Vrs = Vrs[ seq(1, tGap6days*HourScaling, by=HourScaling)];
          
          Vrs_allPeriod = cbind(Vrs_allPeriod, Vrs[1:(24*9)])
          
          Vrs = Vrs[(tGap+1):(tGap6days) ];
          Vrs[which(is.na(Vrs))] = 0
          Vrs_all = cbind(Vrs_all, Vrs);
        }
      }
      
      tempV = rowMeans(Vrs_all);
      
      PlaceboScaleHour12 = c();
      length_min = min(length(tempV), length(temp_PlaceboScaleHour));
      for(i in seq(1,length_min, 24)) {
        PlaceboScaleHour12 = c(PlaceboScaleHour12, mean(temp_PlaceboScaleHour[i:(i+23)] ));
      }
      PlaceboScaleHour12[which(PlaceboScaleHour12<1.1)] = 0.7  
      errorbarVs = c(2.36, 1.96, 1.95);
      
      
      tempV12s = c()
      for (iV in 1:dim(Vrs_all)[2]) {
        tempV = Vrs_all[,iV]
        temp = log10(tempV);
        temp[which(temp<0.7)] = 0.7
        tempV12 = c()
        
        HourMarker = round(runif(1, 0.51, 12.49))
        
        for(iTemp in seq(1, length_min, by=24)) {
          temp24 = 0
          if(iTemp>=24){
            temp24 = round(runif(1, min=1-HourMarker-0.49, max=12-HourMarker+0.49) );
          }
          tempV12= c(tempV12, temp[iTemp + temp24]);
          
        } 
        tempV12s = cbind(tempV12s, tempV12)
      }
      
      
      
      Vrs_all_diff = c()
      for (iV in 1:dim(tempV12s)[2]) {
        Vrs_all_diff = cbind(Vrs_all_diff, tempV12s[,iV]-tempV12s[1,iV])
      }
      
      stdTests_diff= c()
      normalTests = c()
      for (iV in 1:dim(Vrs_all_diff)[1]) {
        tempV = Vrs_all_diff[iV,]
        if( length(unique(tempV))==1 ){
          stdTests_diff = cbind(stdTests_diff, 0)
        }else{
          
          stdTest = sqrt( var(tempV) )
          normalTest = 1
          stdTests_diff = cbind(stdTests_diff, stdTest)
          normalTests = cbind(normalTests, normalTest)
        }
      }
      tempV12 = rowMeans(Vrs_all_diff);
      
      PlaceboScaleHour123 = PlaceboScaleHour12
      PlaceboScaleHour123 = PlaceboScaleHour123 - PlaceboScaleHour123[1]
      
      
      if( length(which(stdTests_diff[c(9)]<1.4)) > 0 || 
          length(which(stdTests_diff[c(2,3,4,5,6,9)]>3)) > 0 ||
          length(which(normalTests[c(2,3,4,5,6,9)]<0.05)) > 0)
      {
        return(-10^10)
      }else{
        if(ei==1){
          y <- dnorm( tempV12, mean = PlaceboScaleHour123, sd = errorbarVs[ei])
          logliks = logliks + sum(log(y))
          
        }else{
          y <- dnorm( tempV12[1:5], mean = PlaceboScaleHour123[1:6], sd = errorbarVs[ei])
          logliks = logliks + sum(log(y))
          
        }
      }
    }
  }
  
  return(logliks)
}




