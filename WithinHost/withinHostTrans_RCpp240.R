## Description:
# Agent-based simulation:
# Simualte each agent, with input of treatment timing since infection, which is the sum of:
# (1) random time for symptom since infection and (2) random time for treatment since symptom
library(Rcpp)

cppFunction('NumericVector withinHostTransRCpp240(int tGap, int tGap6days, double epsilonSet, double betaScaling,
double pScaling, double VsScaling, double cScaling, 
double ZScaling,double  deltaScaling, double rScaling, 
double sigmaScaling, double HourScaling) {
  //double tPeriod = 240;
  double tPeriod = 220;
  NumericVector T(tPeriod+1);
  NumericVector Vs(tPeriod+1);
  NumericVector Vr(tPeriod+1);
  NumericVector Is(tPeriod+1);
  NumericVector Ir(tPeriod+1);
  NumericVector Z(tPeriod+1);
  double dT = 0;
  double dIs = 0;
  double dIr = 0;
  double dZ = 0;
  double dVs= 0;
  double dVr =0;

     Is[0] = 0; 
     Ir[0] = 0; 
     Vr[0] = 0;
     T[0] = 4*pow(10, 8);
     Z[0] = 0.34 * pow(10,ZScaling); //# X0
     Vs[0] = 0.0770 * pow(10,VsScaling); // # Initial sensitive viral load

      double mu = 0;
      double sigma= 0.1*sigmaScaling; //# c
      double delta = 0.0208 * pow(10,deltaScaling); //# dI
      double c = 0.0034 * pow(10,cScaling); // # dV
      
      double beta = 0.0041 * pow(10,betaScaling) ; //# b Infection rate, can adjust the peak timing
      double p = 5 * pow(10,-7) * pow(10,pScaling);// #p Virus production rate, can adjust the peak value
      
      double r = 0.0417 * pow(10,rScaling); //# r Immune response growth rate
      double k = 0.0417;
      double epsilon = 0;

      for (int t = 0; t<tPeriod; t++){
      
         if(t == tGap+1 && epsilonSet>0){
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
            if (Vs[t+1]>Vs[t] && t>48) { Vs[t+1]=Vs[t]; }
            if (Vs[t+1]<0) {Vs[t+1]=0; }
            if (Vr[t+1]<0) {Vr[t+1]=0; }                
                       
      }
      NumericVector Vrs = Vr+Vs;
      return Vrs; 

}')


