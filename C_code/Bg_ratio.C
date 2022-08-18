
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
#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>



std::map<int, std::vector<double>> readtxt(std::string filename){
  std::map<int,std::vector<double> > output;
  std::ifstream infile;
  std::string line;
  infile.open(filename.c_str());
  int i=0;
  while(getline(infile,line)){
    int det;
    double c;  
  
    if(line[0]=='#') continue;  
    if(line[0]==' ') continue;  
    if(line[0]=='\t') continue;  

    std::stringstream ss(line);
    ss >> det;
    while(ss >> c){
      output[det].push_back(c);
    }    
  }

  return output;
}


int main(){

  TFile *f1 = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TFile *f2 = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TFile *fbg = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1115.root");

  std::map<int, std::vector<double>> mat = readtxt("/home/zhu/packages/BCSSort/efficiency/BG_ratio.txt");
  
  TFile *newf = new TFile("efficiency_single_bgsub.root","recreate");
  //TDirectory *current = gDirectory;
  for(int i=0;i<64;i++){
    if(i==7 || i==9 || i==59) continue;
    
    TH1D *h1 = (TH1D *)f1->Get(Form("single_Cal_%i",i));
    TH1D *h2 = (TH1D *)f2->Get(Form("single_Cal_%i",i));
    TH1D *hbg = (TH1D *)fbg->Get(Form("single_Cal_%i",i));
    
    TH1D *ch1 = (TH1D *)h1->Clone(Form("1142_single_%i",i));
    TH1D *ch2 = (TH1D *)h2->Clone(Form("1143_single_%i",i));
    ch1->SetTitle(Form("1142_single_%i",i));    
    ch2->SetTitle(Form("1143_single_%i",i));    

    if(mat[i].size()>2){
      printf("stripe = %i\n",i);
    }
    ch1->Add(hbg,-mat[i][0]); 
    ch2->Add(hbg,-mat[i][1]); 

    //current->cd();
    ch1->Write();
    ch2->Write();
    
  }

  newf->Close();
  return 0;

}
