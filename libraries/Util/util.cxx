


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

//void ReadDetMap(std::string filename) {
void ReadDetMap(const char* filename,Option_t *opt) {

  std::string infilename = filename;

  std::ifstream infile;
  std::string line;
  //infile.open("detmap.txt");
  infile.open(infilename.c_str());
  while(getline(infile,line)) {
    //printf("%s\n",line.c_str());
    int address;
    int crate;
    int slot;
    int chan;
    std::string name;
    int number;
    int segment;
    double offset;
    double gain;
    double quad;
    double time;

    if(line.length()<1) continue;
    if(line[0]=='#') continue;
    if(line[0]==' ') continue;
    if(line[0]=='\t') continue;
    std::stringstream ss(line);
    ss << std::hex;
    ss >> address;

    ss << std::dec;
    ss >> crate;
    ss >> slot;
    ss >> chan;
    ss >> name;
    ss >> number;
    ss >> segment;
    ss >> offset;
    ss >> gain;
    ss >> quad;
    ss >> time;

    TChannel *c = new TChannel(name.c_str()); //TChannel::GetChannel(address);
    //c->SetName(name.c_str());
    c->SetAddress(address);
    c->SetNumber(number);
    c->SetSegment(segment);
    std::vector<double> output;
    output.push_back(offset);
    output.push_back(gain);
    //if(number>=192 || number <=206) {
    output.push_back(quad);
    //}
    c->SetEnergyCoeff(output);

    std::vector<double> time_coeff;
    time_coeff.push_back(time);
    time_coeff.push_back(1.00000);

    c->SetTimeCoeff(time_coeff);

    TChannel::AddChannel(c);
  }
  printf("channels read: %i\n",TChannel::Size());

}


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


void FillHistogram(std::string hname,TCutG *cut,int xbins, double xlow, double xhigh, double xvalue,
                                      int ybins, double ylow, double yhigh, double yvalue) {
    std::string newname = Form("%s_%s",hname.c_str(),cut->GetName());
    FillHistogram(newname,xbins,xlow,xhigh,xvalue, ybins,ylow,yhigh,yvalue);
}

void FillHistogram(std::string hname,int xbins, double xlow, double xhigh, double xvalue,
                                      int ybins, double ylow, double yhigh, double yvalue) {
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
  if(ybins>0) {
    hist->Fill(xvalue,yvalue);
  } else {
    hist->Fill(xvalue);
  }
}


void SaveHistograms(std::string fname) {
  if(!gList) {
    printf("no histograms to save!\n");
    return;
  }
  TDirectory *current = gDirectory;
  TFile *file = new TFile(fname.c_str(),"recreate");
  gList->Sort();
  gList->Write();
  file->Close();
  current->cd();
}



