#ifndef __CINT__
#include <iostream>
#include <numeric>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>


//#include <TApplication.h>
#include <TRint.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>

using namespace std;
#endif


std::map< int,std::vector<double> > readtxt(std::string filename){
 
  std::map<int,std::vector<double> > output;
  ifstream infile;
  string line;
  infile.open(filename.c_str());
  int i=0;
  while(getline(infile,line)){
    int det;
    double c;  
    stringstream ss(line);
    ss >> det;
    while(ss >> c){
      output[det].push_back(c);
    }    
  }

  return output;

}


double GammaEff(double *x, double *par){
  //double logE = TMath::Log(x[0]);
  //double temp =  par[0] + par[1]*logE + par[2]*logE*logE + par[3]*pow(logE,3) +par[4]/(x[0]*x[0]);
  double logE = TMath::Log10(x[0]);
  double temp =  par[0] + par[1]*logE + par[2]*logE*logE + par[3]/(x[0]*x[0]);
  return pow(10,temp);
  //return pow(TMath::E(),temp);
}


TF1 *SRM_Sum(string filename, int det=0){
  std::map<int, std::vector<double>> count = readtxt(filename);
  

  int size = count[0].size();
  double ratio[size];
  double energy[size];
  double err[size] = {2.120788287, 1.313796825, 1.202523222, 1.202523222, 1.202523222, 1.255413119, 1.376975707, 1.255413119, 1.155881524};//err(%) from emitting counts
  double e_err[size];
  
  for(int i=0;i<size;i++){
    energy[i] = count[-1][i]; //count[-1] = peak energy;
    e_err[i] = 0;
    ratio[i] = count[det][i]/count[-2][i]; // count[-2] = emitting counts for each energy; 
    double c_err = 1; // sys_err(%)
    c_err = sqrt(1./count[det][i] + pow((c_err/100.),2)); // calculate error for peak fitting counts;
    err[i] = ratio[i] *sqrt(pow(c_err,2) + pow((err[i]/100.),2)); 
  }
  
  TGraphErrors *gr = new TGraphErrors(size,energy,ratio,e_err,err);
  gr->SetTitle(Form("SRM Sum %i",det));
  TF1 *fxSRM = new TF1("fxSRM",GammaEff,0,4000,4);
  //fxSRM->SetParameters(-1,1,1,-1);
  gr->Fit(fxSRM);
    
  new TCanvas;
  gr->Draw("A*");
  
  for(int i=0;i<size;i++){
    double cal = fxSRM->Eval(energy[i]);
    double dif = fabs((ratio[i]-cal)/cal)*100;
    printf("[%.1f] Observed = %f | Calculated = %f | per dif = %f\n", energy[i],ratio[i],cal,dif);
  }

  return fxSRM;
}


TF1 *Co56_Sum(string filename1, string filename2, int det=0){
  
  TF1 *fx = SRM_Sum(filename1, det);
  std::map<int, std::vector<double>> count = readtxt(filename2);
  double ratio1 = fx->Eval(846.771);
  double ratio2 = fx->Eval(1037.84);
  double ratio3 = fx->Eval(1238.282);
  
  double emit1 = count[det][0]/count[-2][0]/ratio1;
  double emit2 = count[det][2]/count[-2][2]/ratio2;
  double emit3 = count[det][4]/count[-2][4]/ratio3;
  double emit_tot = (emit1 + emit2 + emit3)/3.;

  int size = count[-2].size(); // count[-2] = relative intensity for each peak;
  double ratio[size];
  double err[size];
  double energy[size];
  double e_err[size];
  vector<double> fac;
  
  for(int i=0;i<size;i++){
    double emit = emit_tot * count[-2][i]/100;
    energy[i] = count[-1][i]; //count[-1] = peak energy;
    e_err[i] = 0;
    ratio[i] = count[det][i]/emit;
    if(i<5){
      double factor = fx->Eval(energy[i])/ratio[i];
      fac.push_back(factor);
    } 
    double c_err = 1; // sys_err(%)
    if(c_err < 3){c_err = 3;} // if sys_err < 3% => sys_err = 3%;
    c_err = sqrt(1./count[det][i] + pow((c_err/100.),2)); // calculate error for peak fitting counts;
    err[i] = ratio[i]*c_err;
  } 
  
  double fac_ave = accumulate(fac.begin(), fac.end(),0.0)/fac.size();
  for(int i=0;i<size;i++){
    ratio[i] = ratio[i]*fac_ave;
    err[i]   = err[i]*fac_ave;
  }
  
  TGraphErrors *gr = new TGraphErrors(size,energy,ratio,e_err,err);
  gr->SetTitle(Form("56Co Sum %i",det));
  TF1 *fxCo = new TF1("fxCo",GammaEff,0,4000,5);
  fxCo->SetParameters(ratio[0],0.5,-0.1,2000,1000);
  gr->Fit(fxCo);
    
  new TCanvas;
  gr->Draw("A*");
 
  return fxCo; 

}

void Sum(string filename1, string filename2, int det=0){

  TF1 *f = SRM_Sum(filename1, det);
  std::map<int, std::vector<double>> SRMcount = readtxt(filename1);
  std::map<int, std::vector<double>> SRMsyserr = readtxt("/home/zhu/packages/BCSSort/efficiency/SRM_Sum/SRM_Sum_xtal_syserr.dat");
  std::map<int, std::vector<double>> Cocount = readtxt(filename2);
  std::map<int, std::vector<double>> Cosyserr = readtxt("/home/zhu/packages/BCSSort/efficiency/Co56_Sum/Co56_Sum_xtal_syserr.dat");
  
  double ratio1 = f->Eval(846.771);
  double ratio2 = f->Eval(1037.84);
  double ratio3 = f->Eval(1238.282); 
  double emit1 = Cocount[det][0]/Cocount[-2][0]/ratio1;
  double emit2 = Cocount[det][2]/Cocount[-2][2]/ratio2;
  double emit3 = Cocount[det][4]/Cocount[-2][4]/ratio3;
  double emit_tot = (emit1 + emit2 + emit3)/3.;

  int SRMsize = SRMcount[0].size();
  int Cosize  = Cocount[0].size();
  int size = SRMsize + Cosize;
  double ratio[size];
  double err[size];
  double energy[size];
  double e_err[size]; 
  double c_err = 0;
  vector<double> fac;

  for(int i=0;i<size;i++){
    if(i<SRMsize){
      energy[i] = SRMcount[-1][i];
      e_err[i]  = 0;
      ratio[i]  = SRMcount[det][i]/SRMcount[-2][i];
      c_err = 1; // sys_err(%)
      c_err = sqrt(1./SRMcount[det][i] + pow((c_err/100.),2)); // calculate error for peak fitting counts;
      err[i] = ratio[i] *sqrt(pow(c_err,2) + pow((SRMsyserr[-2][i]/100.),2));  
    }else{
      int j = i-SRMsize;
      double emit = emit_tot * Cocount[-2][j]/100;
      energy[i] = Cocount[-1][j];
      e_err[i]  = 0;
      ratio[i] = Cocount[det][j]/emit;
      if(j<5){
        double factor = f->Eval(energy[i])/ratio[i];
        fac.push_back(factor);
      }
      c_err = 1; // sys_err(%)
      c_err = sqrt(1./Cocount[det][j] + pow((c_err/100.),2)); // calculate error for peak fitting counts;
      err[i] = ratio[i] * c_err;  
    }
  
  }

  double fac_ave = accumulate(fac.begin(), fac.end(), 0.0)/fac.size();
  for(int i=SRMsize;i<size;i++){
    ratio[i] = ratio[i] * fac_ave;
    err[i]   = err[i] *fac_ave; 
  }  

  //for(int i=0;i<size;i++) cout<<energy[i] <<"\t" <<ratio[i] <<"\t" << err[i]<<endl;

  TGraphErrors *gr = new TGraphErrors(size,energy,ratio,e_err,err);
  gr->SetTitle(Form("Sum %i", det));
  //gr->SetTitle(Form("Area %i", det));
  TF1 *fx = new TF1("fx",GammaEff,0,4000,5);
  fx->SetParameters(-2,0.3,-0.2,0.1,1000);
  gr->Fit(fx);
    
  new TCanvas;
  gr->Draw("A*");

  for(int i=0;i<size;i++){
    double cal = fx->Eval(energy[i]);
    double dif = fabs((ratio[i]-cal)/cal)*100;
    printf("[%.1f] Observed = %f | Calculated = %f | per dif = %f\n", energy[i],ratio[i],cal,dif);
  }

  return;
}


int main(int argc, char **argv){
  TRint app("app",0,0);
  std::string arg;
  for(int i=1;i<argc;i++){
    arg = argv[i];
    //if(isdigit(arg)==false){
      SRM_Sum(arg,0);
    //} 
  }
  app.Run(0);
  
  return 0;
}

/*int main(int argc, char **argv){
  //TApplication app("app",0,0);
  TRint app("app",0,0);
  
  int det = 0;
  string arg;
  if(argc<3){ 
    det = 0;
    arg = argv[1];
    Sum(argv[1],det);
  }
  else{
    int i=1;
    arg = argv[2];
    if(isdigit(arg)) {
      for(i=3;i<argc;i++){
        arg = argv[i];
        det = std::stoi(arg);
        Sum(argv[1],argv[2],det); 
      }
    }
    else {
      for(i=2;i<argc;i++){
        arg = argv[i];
        det = std::stoi(arg);
        SRM_Sum(argv[1],det);
      }
    }
  }

  app.Run(0);

  return 0;
}*/


















