
#include <OutputManager.h>

#include <iostream>
#include <cstdio>

#include <TFile.h>
#include <TTree.h>

#include <Implant.h>
#include <BCSOptions.h>

OutputManager *OutputManager::fOutputManager =0;


OutputManager::OutputManager() : fFile(0), fTree(0), fBeta(0), fEventFile(0), fEventTree(0), fEvent(0), fListFile(0) ,fListTree(0), faddress(0), fnumber(0), ftimestamp(0), fcharge(0){ } 

OutputManager::~OutputManager() { }

OutputManager *OutputManager::Get() {
  if(!fOutputManager) {
    fOutputManager = new OutputManager;
  }
  return fOutputManager;
}


const char *OutputManager::GetName() const { if(fFile) { return fFile->GetName(); } return ""; }

void OutputManager::Set(std::string runnumber) {
  if(BCSOptions::Get()->WriteTree()) {
   
    TDirectory * current = gDirectory; 
    
    fFile = new TFile(Form("beta%s.root",runnumber.c_str()),"recreate");
    fListTree = new TTree("beta","beta tree");
    fListTree->SetDirectory(fFile);
    fBeta = new Beta;
    fListTree->Branch("Beta",&fBeta);  
   
    current->cd();
  }
  if(BCSOptions::Get()->WriteEventTree()) {
    TDirectory * current = gDirectory; 
    
    fEventFile = new TFile(Form("event%s.root",runnumber.c_str()),"recreate");
    fEventTree = new TTree("event","event tree");
    fEventTree->SetDirectory(fEventFile);
    fEvent  = new BCSEvent;
    fEventTree->Branch("BCSEvent",&fEvent);  
   
    current->cd();

  }
  if(BCSOptions::Get()->WriteListTree()) {
    TDirectory * current = gDirectory; 
    
    fListFile = new TFile(Form("list%s.root",runnumber.c_str()),"recreate");
    fListTree = new TTree("tree","tree");
    fListTree->SetDirectory(fListFile);
    //fHit  = new DetHit;
    //fListTree->Branch("DetHit",&fHit);
    //faddress = new int;
    //fnumber = new int;
    //ftimestamp = new double;
    //fcharge = new double;
    fListTree->Branch("address",    &faddress,    "address/I");  
    fListTree->Branch("number",     &fnumber,     "number/I");  
    fListTree->Branch("timestamp",  &ftimestamp,  "timestamp/D");  
    fListTree->Branch("charge",     &fcharge,     "energy/D");  
   
    current->cd();

  }
  return;
}


void OutputManager::Fill(Implant *implant,std::vector<Decay> *vdec) {
  if(!BCSOptions::Get()->WriteTree()) return;
  
  //TDirectory * current = gDirectory; 
  //fFile->cd();
  fBeta->Set(*implant, *vdec);
  //fBeta->SimplePrint();
  fListTree->Fill();
  fBeta->Clear();

  //current->cd();
}


void OutputManager::FillEvent(std::vector<DetHit> *hits) {
  if(!BCSOptions::Get()->WriteEventTree()) return;

  fEvent->fHits = *hits;
  fEventTree->Fill();
  fEvent->Clear();

}

void OutputManager::FillList(int address, int number, double timestamp, double charge) {
  //if(!BCSOptions::Get()->WriteListTree()) return;

   faddress   = address;
   fnumber    = number;
   ftimestamp = timestamp;
   fcharge    = charge;
   fListTree->Fill();

}

void OutputManager::Close() {
  //if(!BCSOptions::Get()->WriteTree()) return;
  if(BCSOptions::Get()->WriteTree()) {
    if(!fFile || !fListTree) return;
    TDirectory * current = gDirectory; 
    fFile->cd();
    fListTree->Write();
    fFile->Close();
    current->cd();
  }
  if(BCSOptions::Get()->WriteEventTree()) {
    if(!fEventFile || !fEventTree) return;
    TDirectory * current = gDirectory; 
    fEventFile->cd();
    fEventTree->Write();
    fEventFile->Close();
    current->cd();
  }
  if(BCSOptions::Get()->WriteListTree()) {
    if(!fListFile || !fListTree) return;
    TDirectory * current = gDirectory; 
    fListFile->cd();
    fListTree->Write();
    fListFile->Close();
    current->cd();
  }


  return;
}

