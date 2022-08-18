

#include<cstdio>
#include<TGraph.h>
#include<TGraphErrors.h>


//================================= Fitting Model ==================================//

// simple exponential function;
// par[0]: constant scaling factor
// par[1]: tau (halflife)
// dim[0]: x value - time 
double MYExp(double *x, double *par){

  double lam1 = 0.693/par[1];
  return par[0]*TMath::Exp(-lam1*x[0]);

}

// Calculate activity of one decay = lambda *N(t);
// N(t) is calculated by Bateman equation;
// n = isotope number(parent = 1, daughter = 2, granddaughter = 3);
// when n = 1, return exponential function par[0]*exp(-lam1*x[0])
// par[0] = particle[n] activity;par[i] = half-life of particle_i(i>0)
double Bateman(double *x, double *par, int n=1){
  double mul = 1;
  double sum = 0;
  for(int i=1;i<=n;i++){
    double lami = 0.693/par[i];
    if(i>=2 && (i<=n-1)) mul = mul*lami;
    double denominator = 1;
    for(int j=1;j<=n;j++){
      if(j==i) continue;
      double lamj = 0.693/par[j];
      denominator = denominator * (lamj - lami);
    }
    sum += TMath::Exp(-lami*x[0])/denominator; 
  }
  if(n==1) return par[0]*mul*sum;
  return (0.693/par[n])*(par[0]*mul*sum);
}



// sum of activities for a decay chain;
// n =  # of isotopes decay in this decay chain;
// # of par = n+1(n half-lifes + 1 parent activity);
// par[0] = the only one activity for each chain
double sum_Bateman(double *x, double *par, int n){
  double sum = 0;
  for(int i=1;i<=n;i++){
    sum += Bateman(x,par,i);
  }
  return sum;
}

//Eventual Fitting Model
//Example: 2 decay paths: 32Na->32Mg->32Al & 32Na->31Mg
//two bg: constant + exp bg
double MYFit(double *x, double *par) {
  //par[0]: constant bg;
  //par[1]: A_bg
  //par[2]: t1/2_bg;
  //par[3]: A: path1: 30Ne->30Na->30Mg->30Al;
  //par[4]: A: path2: 30Ne->30Na->29Mg;
  //par[5]: A: path3: 30Ne->29Na->29Mg
  //par[6]: t1/2: 30Ne;
  //par[7]: t1/2: 30Na;   
  //par[8]: t1/2: 29Na;
  //par[9]: t1/2: 30Mg;

 
  double par_bg[2];
  double par1[4]; //path1: 30Ne->30Na->30Mg->30Al
  double par2[3]; //path2: 30Ne->30Na->29Mg
  double par3[3]; //path3: 30Ne->29Na->29Mg

  par_bg[0] = par[1];
  par_bg[1] = par[2];

  par1[0] = par[3];
  par1[1] = par[6];
  par1[2] = par[7];
  par1[3] = par[9];

  par2[0] = par[4];
  par2[1] = par[6];
  par2[2] = par[7];
  
  par3[0] = par[5];
  par3[1] = par[6];
  par3[2] = par[8];
  
  double value = par[0] + MYExp(x,par_bg)+ sum_Bateman(x,par1,3) + sum_Bateman(x,par2,2) + sum_Bateman(x,par3,2);
  
  return value;
}
//================================== Math model for fitting plot =================================//
//par[0] = 0 => bg plot;
//par[0] >=0 => the value of par[0] = isotope number (parent=1, daughter=2, granddaughter=3)
double Draw_fit(double *x, double *par){
  double sum = 0;
  if(par[0]==0) sum = par[1] + MYExp(x,par+2);
  else sum = par[1] + Bateman(x,par+2,par[0]);
  return sum;
}



