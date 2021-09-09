
#include <OutputManager.h>

#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include <Implant.h>
#include <BCSOptions.h>

OutputManager *OutputManager::fOutputManager =0;


OutputManager::OutputManager() : fFile(0), fTree(0), fBeta(0), fEventFile(0), fEventTree(0), fEvent(0), fListFile(0) ,fHitTree(0), fHit(0){ } 

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
    fTree = new TTree("beta","beta tree");
    fTree->SetDirectory(fFile);
    fBeta = new Beta;
    fTree->Branch("Beta",&fBeta);  
   
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
    fHitTree = new TTree("hit_tree","hit tree");
    fHitTree->SetDirectory(fListFile);
    fHit  = new DetHit;
    fHitTree->Branch("DetHit",&fHit);  
   
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
  fTree->Fill();
  fBeta->Clear();

  //current->cd();
}


void OutputManager::FillEvent(std::vector<DetHit> *hits) {
  if(!BCSOptions::Get()->WriteEventTree()) return;

  fEvent->fHits = *hits;
  fEventTree->Fill();
  fEvent->Clear();

}

void OutputManager::FillList(DetHit *hit) {
  //if(!BCSOptions::Get()->WriteListTree()) return;

  fHit = hit;
  fHitTree->Fill();

}

void OutputManager::Close() {
  //if(!BCSOptions::Get()->WriteTree()) return;
  if(BCSOptions::Get()->WriteTree()) {
    if(!fFile || !fTree) return;
    TDirectory * current = gDirectory; 
    fFile->cd();
    fTree->Write();
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
    if(!fListFile || !fHitTree) return;
    TDirectory * current = gDirectory; 
    fListFile->cd();
    fHitTree->Write();
    fListFile->Close();
    current->cd();
  }


  return;
}

