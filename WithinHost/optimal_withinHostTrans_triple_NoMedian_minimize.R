## Description:
# Agent-based simulations
# (1) Simualte each agent, with (1) random time for symptom since infection and (2) random time for treatment since symptom
# (2) Estimate MCMC likelihood as Joint Normal probability of agents with observed mean and stds.

optimal_withinHostTrans_triple_NoMedian_minimize= function(x, PopDist_Baloxavir, BaloxavirScaleHour,HourScaling,
                                                           PopDist_Placebo, PlaceboScaleHour, PopDist_Oseltamivir, 
                                                           OseltamivirScaleHour, lb, ub, returnLog, eiSingle, parasNum, otherParas){
  
  logliks =  optimal_withinHostTrans_triple_Median_Symptom(x, PopDist_Baloxavir, BaloxavirScaleHour,HourScaling,
                                                           PopDist_Placebo, PlaceboScaleHour, PopDist_Oseltamivir, 
                                                           OseltamivirScaleHour, lb, ub, returnLog, eiSingle, parasNum, otherParas)
  return(-logliks)
}




