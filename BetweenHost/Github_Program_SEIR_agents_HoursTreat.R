
Program_SEIR_agents_HoursTreat = function(beta,sigma, gamma, S0, E0, I0, N0, MaxTime, itry, pVisit, drugEffiV_Bt, pDayOne, xOthers){
  
  set.seed(itry)
  
  S0 = ceil(S0);
  E0 = ceil(E0);
  I0 = ceil(I0);
  R0 = N0-S0-E0-I0;
  D = (1/gamma);
  incidence = Stoch_Iteration(c(0, MaxTime), c(S0, E0, I0, R0), c(beta, gamma, sigma, N0, D), pVisit, drugEffiV_Bt, pDayOne, xOthers);
  
}

Stoch_Iteration = function(Time, Initial, Parameters, pVisit, drugEffiV_Bt, pDayOne, xOthers){
  
  X=Initial[1];
  Y = Initial[3];
  Z = Initial[4];
  T=1;
  # P[1,]=c(X, Y, Z)
  P=c(X, Y, Z)
  old=c(X, Y, Z)
  
  # #%#%
  popSize = Parameters[4];
  List_S = rep(1,popSize)>0; 
  List_E = rep(1,popSize)<0; 
  List_I = rep(1,popSize)<0; 
  List_R = rep(1,popSize)<0; 
  List_B = rep(1,popSize)<0; 
  
  # #%#%
  List_betaIt = rep(0, popSize);
  List_betaBt = rep(0, popSize);
  
  # #% lists of Timings
  List_tI = rep(0, popSize);
  List_tB = rep(0, popSize);
  List_tBhour = rep(0, popSize);
  
  List_pIt = rep(0, popSize);
  List_pBt = rep(0, popSize);
  List_pTreat = rep(0, popSize);
  
  PopDistNum_Hours_Len = c(0,0,0)
  while(length(which(PopDistNum_Hours_Len[5:length(PopDistNum_Hours_Len)]==0))){
    if(pDayOne>0){
      PopDistNum_Hours2= c();
      temp10000 =  round(rgamma(10000, shape=3.9855, rate = 1/6.3190))
      temp24 = which(temp10000<=24)
      temp48 = which(temp10000>24&temp10000<=48)
      PopDistNum_Hours2 = c(temp10000[temp24[1:(pDayOne*1000)]],  temp10000[temp48[1:((1-pDayOne)*1000) ]])
      
      PopDistNum_Hours2 = PopDistNum_Hours2 + round(rlnorm(length(PopDistNum_Hours2), meanlog=0.3686, sdlog = 0.4052)*24)
      
    }
    if(pDayOne==-1){
      # PopDistNum_Hours2= c();
      temp10000 =  round(rgamma(10000, shape=3.9855, rate = 1/6.3190))
      temp24 = which(temp10000<=24)
      temp48 = which(temp10000>24&temp10000<=48)
      PopDistNum_Hours2 = c(temp10000[temp24[1:(0.53*1000)]],  temp10000[temp48[1:((1-0.53)*1000) ]])
      # PopDistNum_Hours2 = round( c(runif(135, 1, 12)-1, runif(428, 1, 12)+12-1, runif(300, 1, 12)-1+24, runif(201, 1, 12)+12-1+24)  )
      PopDistNum_Hours2 = PopDistNum_Hours2 + round(rlnorm(length(PopDistNum_Hours2), meanlog=0.3686, sdlog = 0.4052)*24)
    }
    if(pDayOne<=-2){ # pDayOne = -12, denoting 12 hours treatment all agents
      PopDistNum_Hours2= c();
      PopDistNum_Hours2 = -pDayOne; #12
      PopDistNum_Hours2 = PopDistNum_Hours2 + round(rlnorm(1000, meanlog=0.3686, sdlog = 0.4052)*24)
    }
    
    PopDistNum_Hours = list(c());
    PopDistNum_Hours_Len = c()
    for(dayi in 1:ceil(max(PopDistNum_Hours2)/24)){
      temp = PopDistNum_Hours2[which(PopDistNum_Hours2<24*dayi & PopDistNum_Hours2>=24*(dayi-1) )]
      PopDistNum_Hours[[dayi]] = temp - 24*(dayi-1)
      PopDistNum_Hours_Len = c(PopDistNum_Hours_Len, length(PopDistNum_Hours[[dayi]]) )
    }
  }
  # print(PopDistNum_Hours_Len)
  
  ## %#%
  if(pDayOne>=0){
    Results = getpITreatAdjust(Parameters, pVisit, pDayOne, PopDistNum_Hours_Len);# #%pBt
    pIt= Results[[1]]
    pTreat= Results[[2]]
    pBt = Results[[3]]
  }
  if(pDayOne==-1){
    Results = getpITreat(Parameters, pVisit, pDayOne, PopDistNum_Hours_Len); ##%pBt
    pIt = Results[[1]]
    pTreat= Results[[2]]
    pBt = Results[[3]]
  }
  if(pDayOne<=-2){ # if less than -2, it denots the hours at.
    Results = getpITreatAdjust(Parameters, pVisit, 0.53, PopDistNum_Hours_Len);# #%pBt
    pIt= Results[[1]]
    pTreat= Results[[2]]
    pBt = Results[[3]]
  }
  
  Results= getbetaIBtHours(drugEffiV_Bt, xOthers)
  betaIt = Results[[1]]
  betaIt[is.na(betaIt)]=0
  
  betaBtHours = Results[[2]]
  
  # Results= getbetaIBt(drugEffiV_Bt, xOthers)
  # betaIt = Results[[1]]
  # betaBt = Results[[2]]
  
  # #%#% ini
  library("pracma")
  tempS = which(List_S==1);
  tempI = tempS[randperm(length(tempS), Y)];
  List_I[tempI] = 1
  List_S[tempI] = 0;
  
  List_tI[tempI] = 1;
  List_pIt[tempI] = pIt[1];
  List_betaIt[tempI] = betaIt[1];
  
  tempS2 = which(List_S==1);
  tempR = tempS2[randperm(length(tempS2), Z)];
  List_R[tempR] = 1; List_S[tempR] = 0;
  
  # #%#%
  loop=1
  incidence=rep(0,(Time[2]+1)*2)
  # incidence=zeros((Time[2]+1),1)
  tempTrack=c()
  
  while (T[loop]<Time[2]){
    tNow = T[loop];
    if (sum(List_I)==0){  break; }
    
    Results = Iterate(old, Parameters, List_R, List_B, List_I, List_S, List_E,
                      pIt, pTreat,
                      betaIt, betaBtHours,
                      List_betaIt, List_betaBt,
                      List_tI, List_tB,
                      List_pBt, List_pTreat, List_pIt, tNow, List_tBhour, PopDistNum_Hours);
    step=Results[[1]]
    new=Results[[2]]
    m=Results[[3]]
    List_R=Results[[4]]
    List_B=Results[[5]]
    List_I=Results[[6]]
    List_S=Results[[7]]
    List_E=Results[[8]]
    List_tB=Results[[9]]
    List_tI=Results[[10]]
    List_pBt=Results[[11]]
    List_pIt=Results[[12]]
    List_pTreat=Results[[13]]
    List_betaIt=Results[[14]]
    List_betaBt=Results[[15]]
    
    List_tBhour=Results[[16]]
    
    if(is.infinite(step) || m==-1){
      break
    }
    temp_t = tNow+step;
    Results = Updates(temp_t, pIt, pBt, List_tI, List_I, List_pIt, List_B, List_pBt,
                      List_betaIt, List_betaBt, betaIt, betaBtHours, List_tB, List_pTreat, pTreat, List_tBhour);
    List_pTreat=Results[[1]]
    List_pIt=Results[[2]]
    List_pBt=Results[[3]]
    List_betaIt=Results[[4]]
    List_betaBt=Results[[5]]
    
    loop=loop+1;
    T[loop] = T[loop-1]+step;
    # P[loop,]=old;
    loop=loop+1;
    T[loop] = T[loop-1];
    # P[loop,]=new
    old=new;
    
    if (ceil(T[loop])>length(incidence) ){
      break;
    }
    incidence[ceil(T[loop])] = incidence[ceil(T[loop])]+m;
    if (loop>=length(T) ){
      T[loop*2]=0;
      # P[loop*2,]=0;
    }
  }
  
  # T=T(1:loop);
  # P=P[1:loop,]
  return(incidence)
}

#% Do the actual iteration step
Iterate = function (old, Parameters,
                    List_R, List_B, List_I, List_S, List_E,
                    pIt, pTreat,
                    betaIt, betaBtHours,
                    List_betaIt, List_betaBt,
                    List_tI, List_tB,
                    List_pBt, List_pTreat, List_pIt, tNow, List_tBhour, PopDistNum_Hours)
  #% beta=Parameters[1]; #%gamma=Parameters(2); mu=Parameters(3);
{
  step = 1; new = 0; m = 0;
  phi = Parameters[1];
  gamma = Parameters[2]; #% I->R
  sigma = Parameters[3]; #% E->I
  popSize=Parameters[4];
  
  tempBsum = sum(List_betaBt[which(List_B==1)]);
  List_betaIBt_Sum = sum(List_betaIt[which(List_I==1)])+tempBsum;
  lambda = phi*List_betaIBt_Sum;
  
  for (i in 1:length(List_S) ){
    flag = 0;
    if (List_S[i]==1 && flag == 0){
      temp = lambda/popSize;
      if (is.na(temp)){
        temp = 0
      }
      if (runif(1)<=temp){
        flag = 1;
        List_S[i]=0; List_E[i]=1;
      }
    }
    #%
    if (List_E[i]==1 && flag == 0){
      if (runif(1)<=sigma){
        flag = 1;
        List_E[i]=0; List_I[i]=1;
        List_tI[i] = ceil(tNow);
        # List_pTreat[i] = pTreat[1];
        List_pTreat[i] = 0;
        m = m+1;
      }
    }
    #%
    if (List_I[i]==1 && flag == 0){
      if (runif(1)<= gamma){
        flag = 1;
        List_I[i]=0;
        List_R[i]=1;
        List_pIt[i] = 0;
        List_pTreat[i] = 0;
        # List_pTreat[i] = pTreat[1];
        List_betaIt[i] = 0;
        List_betaBt[i] = 0;
        # m = m+1;
      }
    }
    #%
    if (List_I[i]==1 && flag == 0 && List_pTreat[i]>0){
      temp = List_pTreat[i];
      if (runif(1)<=temp){
        flag = 1;
        List_I[i] = 0;
        List_B[i] = 1;
        List_tB[i] = ceil(tNow);
        
        temp_PopDistNum_Hours = PopDistNum_Hours[List_tB[i]-List_tI[i]]
        List_tBhour[i] = temp_PopDistNum_Hours[[1]][randi(length(temp_PopDistNum_Hours[[1]]))]
        
        List_pIt[i] = 0;
        List_pTreat[i] = 0;
        List_betaIt[i] = 0;
      }
    }
    #%
    if (List_B[i]==1&& flag == 0 && List_pBt[i]>0){
      temp = List_pBt[i];
      if (runif(1)<=temp){
        flag = 1;
        List_B[i]=0;
        List_R[i]=1;
        List_pBt[i] = 0;
        List_pIt[i] = 0;
        List_pTreat[i] = 0;
        List_betaIt[i] = 0;
        List_betaBt[i] = 0;
      }
    }
  }
  return(list( step, new, m, List_R, List_B, List_I, List_S, List_E,
               List_tB, List_tI, List_pBt, List_pIt,  List_pTreat,
               List_betaIt, List_betaBt, List_tBhour ))
}

Updates = function(temp_t, pIt, pBt, List_tI, List_I, List_pIt, List_B, List_pBt, List_betaIt, List_betaBt,
                   betaIt, betaBtHours, List_tB, List_pTreat, pTreat, List_tBhour)
{
  temp_List_I = which(List_I==1);
  if(length(temp_List_I)>0){
    for(i in 1:length(temp_List_I) ){
      temp_i = temp_List_I[i];
      if (ceil(temp_t) != ceil( List_tI[temp_i] ) ){
        if (ceil(temp_t)-ceil( List_tI[temp_i] )>0){
          if (ceil(temp_t)-ceil( List_tI[temp_i] )<=length(pIt)){
            List_pIt[temp_i] = pIt[ceil(temp_t)-ceil( List_tI[temp_i] )  ];
          }else{
            List_pIt[temp_i] = 0;
          }
        }
        else{
          List_pIt[temp_i] = 0;
        }
      }
      
      if (ceil(temp_t)-ceil( List_tI[temp_i] )+1 <= length(betaIt) ){
        List_betaIt[temp_i] = betaIt[ ceil(temp_t)-ceil( List_tI[temp_i] )+1 ];}
      else{
        List_betaIt[temp_i] = 0;
      }
      
      if (ceil(temp_t)-ceil( List_tI[temp_i] )+1 <= length(pTreat)){
        List_pTreat[temp_i] = pTreat[ ceil(temp_t)-ceil( List_tI[temp_i] ) ];#%+1 );
      }else{
        List_pTreat[temp_i] = 0;
      }
    }
  }
  
  temp_List_B = which(List_B==1);
  if(length(temp_List_B)>0){
    for( i in 1:length(temp_List_B) ){
      temp_b = temp_List_B[i];
      if (ceil(temp_t)-ceil( List_tI[temp_b] )+1 <= dim(betaBtHours)[1] &&
          ceil(temp_t)-ceil( List_tB[temp_b] )+1 <= dim(betaBtHours)[3]) {
        # temp_Hours = round(rlnorm(1, meanlog=0.3686, sdlog = 0.4052)*24)  +  List_tBhour[temp_b]+1
        # betaBtHours1 = floor(temp_Hours/24)
        # betaBtHours3 = temp_Hours-betaBtHours1*24+1
        
        
        # print(c(ceil(temp_t)-ceil( List_tB[temp_b] )+betaBtHours1,
        #         ceil(temp_t)-ceil( List_tI[temp_b] ), betaBtHours3))
        # List_betaBt[temp_b] = betaBtHours[ ceil(temp_t)-ceil( List_tB[temp_b] ),
        #                                    List_tBhour[temp_b]+1, ceil(temp_t)-ceil(List_tI[temp_b]) ]
        List_betaBt[temp_b] = betaBtHours[ ceil(temp_t)-ceil(List_tI[temp_b]) ,
                                           List_tBhour[temp_b]+1, ceil(temp_t)-ceil( List_tB[temp_b] )]
        # List_betaBt(temp_b) = betaBt( ceil(temp_t)-ceil( List_tI(temp_b) ),
        #                               ceil(temp_t)-ceil( List_tB(temp_b) ), List_tBhour(temp_b)+1 );
      }
      else{
        List_betaBt[temp_b] = 0;
      }
      
      List_pTreat[temp_b] = 0;
      List_betaIt[temp_b] = 0;
      List_pIt[temp_b] = 0;
      #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
      if (ceil(temp_t) != ceil( List_tB[temp_b] ) ){
        if (ceil(temp_t)-ceil( List_tB[temp_b] )>0 ){
          if (ceil(temp_t)-ceil( List_tB[temp_b] )<=length(pBt)){
            List_pBt[temp_b] = pBt[ceil(temp_t)-ceil( List_tB[temp_b] )];
          }else{
            List_pBt[temp_b] = 0;
          }
        }else{
          List_pBt[temp_b] = 0;
        }
      }else{
        List_pBt[temp_b] = pBt[1];
      }
    }
  }
  return(list(List_pTreat, List_pIt, List_pBt, List_betaIt, List_betaBt))
}