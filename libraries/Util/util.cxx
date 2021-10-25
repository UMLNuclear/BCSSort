


#include "util.h"

#include <fstream>
#include <sstream>

#include "TChannel.h"

#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPRegexp.h>

#include <TDirectory.h>
#include <TList.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


TList *gList = new TList;



std::string GetRunNumber(std::string input) {
  TPRegexp re(
      "(" // Begin capturing group.  This will hold the final return value
      "([0-9]+(-|_))?" //One or more digits, followed by a dash or an underscore.  This section may be omitted.
      "[0-9]+" //Followed by one or more digits
      ")" // End capturing group.
      "[^0-9]*$" // With no other digits before the end of the filename
      );
  TObjArray* matches = re.MatchS(input.c_str());

  std::string output;
  //std::cout << " matches->GetEntriesFast() = " << matches->GetEntriesFast() << std::endl;
  //for(int x=0;x<matches->GetEntriesFast();x++)
  //  std::cout << x << "\t" << ((TObjString*)matches->At(x))->GetString() << std::endl;

  if(matches->GetEntriesFast() >= 2){
    // Return a std::vector<std::string> ?
    // No, that would be too simple, too type-safe, and too memory-safe for ROOT.
    output = ((TObjString*)matches->At(1))->GetString();
  }
  delete matches;

  return output;
}


void ProgressUpdate(long x, long n, bool end) {
  if(!end) {
    if((x%50000)==0) {
      printf(" sorting entry %lu / %lu        \r",x,n);
      fflush(stdout);
    }
  } else {
    printf(" sorting entry %lu / %lu        \n",x,n);
  }
}


void FillHistogram(std::string hname,TCutG *cut1, TCutG *cut2, int xbins, double xlow, double xhigh, double xvalue,
                                                                int ybins, double ylow, double yhigh, double yvalue, 
                                                                double zvalue) {
  std::string newname = Form("%s_%s_%s",hname.c_str(),cut1->GetName(),cut2->GetName());
  FillHistogram(newname,xbins,xlow,xhigh,xvalue, ybins,ylow,yhigh,yvalue,zvalue);
}

void FillHistogram(std::string hname,TCutG *cut,int xbins, double xlow, double xhigh, double xvalue,
                                                int ybins, double ylow, double yhigh, double yvalue, 
                                                double zvalue) {
  std::string newname = Form("%s_%s",hname.c_str(),cut->GetName());
  FillHistogram(newname,xbins,xlow,xhigh,xvalue, ybins,ylow,yhigh,yvalue,zvalue);
}


void FillHistogram(std::string hname,int xbins, double xlow, double xhigh, double xvalue,
                                      int ybins, double ylow, double yhigh, double yvalue, 
                                      double zvalue) {

  if(!gList) gList = new TList;

  TH1 *hist = (TH1*)gList->FindObject(hname.c_str());


  if(!hist) {

    if(ybins>0) {
      hist = new TH2D(hname.c_str(),hname.c_str(),xbins,xlow,xhigh,ybins,ylow,yhigh);
    } else {
      hist = new TH1D(hname.c_str(),hname.c_str(),xbins,xlow,xhigh);
    }
    gList->Add(hist);
  }
  
  if(ybins>0 && zvalue==zvalue) {
    ((TH2*)hist)->Fill(xvalue,yvalue,zvalue); // fill with "weighted" zvalue.
    //hist->SetBinContent(xbin,ybin,zvalue);
  } else if(ybins>0) {
    hist->Fill(xvalue,yvalue);
  } else {
    hist->Fill(xvalue);
  }
}

void SaveHistograms(std::string fname, Option_t *opt) {
  if(!gList) {
    printf("no histograms to save!\n");
    return;
  }
  TString sopt(opt);
  sopt.ToLower();
  TDirectory *current = gDirectory;
  TFile *file;
  if(sopt.Length()==0 || sopt.Contains("recreate")){
    file = new TFile(fname.c_str(),"recreate");
  }
  if(sopt.Contains("update")){
    file = new TFile(fname.c_str(),"update");
  }
  gList->Sort();
  gList->Write();
  file->Close();
  current->cd();
}



