
#include "BCSintFunctions.h"

#include <cstdio>

#include <TROOT.h>
#include <TFile.h>

#include <BCSint.h>

void SendHelp() {
   printf("PLEASE SEND HELP!\n");
 }

void OpenRootFile(std::string fname) {

  std::string rfile = Form("_file%i",BCSint::Get()->UpdateFileCount());
  //gROOT->ProcessLine(Form("TFile *%s = new TFile(\"%s\");",rfile.c_str(),fname.c_str()));  
  BCSint::Get()->ProcessLine(Form("TFile *%s = TFile::Open(\"%s\");",rfile.c_str(),fname.c_str()));
  std::cout << "file " << fname << " opened as " << rfile << std::endl;

  //TFile *f = (TFile*)gROOT->FindObjectAny(rfile.c_str());
  TFile *f= (TFile*)gROOT->GetListOfFiles()->Last();
  //BCSint::Get()->ProcessLine(Form("f = %s;",rfile.c_str()));
  //std::cout << f->GetName() << std::endl;
  //if(f->FindObjectAny("dchan")) {
  //if(f->FindObjectAny("tree")) {
  if(f->FindObjectAny("event")) {
  //if(f->FindObjectAny("beta")) {
    gChain->Add(fname.c_str());
  }

}

void LoadCuts() {
    TFile *mycut = TFile::Open("~zhu/notebooks/Ne31/Time_Correlation/mynewcuts.root");
    TIter keys (mycut->GetListOfKeys());
    while(TKey *key = (TKey*)keys.Next()) {
        gCuts->Add(key->ReadObj());
    }

}


