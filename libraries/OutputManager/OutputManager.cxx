
#include <OutputManager.h>

#include <iostream>
#include <cstdio>

#include <TFile.h>
#include <TTree.h>

#include <Implant.h>
#include <BCSOptions.h>

OutputManager *OutputManager::fOutputManager =0;


OutputManager::OutputManager() : fFile(0), fTree(0), fBeta(0), 
                                  fEventFile(0), fEventTree(0), fEvent(0), 
                                  fListFile(0) ,fListTree(0), faddress(0), fnumber(0), ftimestamp(0), fcharge(0),
                                  fImpFile(0), fImpTree(0), fImp(0), 
                                  fDecFile(0), fDecTree(0), fDec(0){ } 

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
    fListTree = new TTree("tree","tree");
    fListTree->SetDirectory(fListFile);
    fListTree->Branch("address",    &faddress,    "address/I");  
    fListTree->Branch("number",     &fnumber,     "number/I");  
    fListTree->Branch("timestamp",  &ftimestamp,  "timestamp/D");  
    fListTree->Branch("charge",     &fcharge,     "energy/D");  
   
    current->cd();

  }
  if(BCSOptions::Get()->WriteImpTree()) {
    TDirectory * current = gDirectory; 
    
    fImpFile = new TFile(Form("implant%s.root",runnumber.c_str()),"recreate");
    fImpTree = new TTree("implant","implant tree");
    fImpTree->SetDirectory(fImpFile);
    fImp  = new Implant;
    fImpTree->Branch("Implant",&fImp);  
   
    current->cd();

  }
  if(BCSOptions::Get()->WriteDecTree()) {
    TDirectory * current = gDirectory; 
    
    fDecFile = new TFile(Form("decay%s.root",runnumber.c_str()),"recreate");
    fDecTree = new TTree("decay","decay tree");
    fDecTree->SetDirectory(fDecFile);
    fDec  = new Decay;
    fDecTree->Branch("Decay",&fDec);  
   
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

void OutputManager::FillList(int address, int number, double timestamp, double charge) {
  //if(!BCSOptions::Get()->WriteListTree()) return;

   faddress   = address;
   fnumber    = number;
   ftimestamp = timestamp;
   fcharge    = charge;
   fListTree->Fill();

}

void OutputManager::FillImp(std::vector<DetHit> *hits){
  if(!BCSOptions::Get()->WriteImpTree()) return;
  
  fImp->Set(hits);
  fImpTree->Fill();
  fImp->Clear();
}
void OutputManager::FillDec(std::vector<DetHit> *hits){
  if(!BCSOptions::Get()->WriteDecTree()) return;
  
  fDec->Set(hits);
  fDecTree->Fill();
  fDec->Clear();
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
    if(!fListFile || !fListTree) return;
    TDirectory * current = gDirectory; 
    fListFile->cd();
    fListTree->Write();
    fListFile->Close();
    current->cd();
  }
  if(BCSOptions::Get()->WriteImpTree()) {
    if(!fImpFile || !fImpTree) return;
    TDirectory * current = gDirectory; 
    fImpFile->cd();
    fImpTree->Write();
    fImpFile->Close();
    current->cd();
  }
  if(BCSOptions::Get()->WriteDecTree()) {
    if(!fDecFile || !fDecTree) return;
    TDirectory * current = gDirectory; 
    fDecFile->cd();
    fDecTree->Write();
    fDecFile->Close();
    current->cd();
  }


  return;
}

