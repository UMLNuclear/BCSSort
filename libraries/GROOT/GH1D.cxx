
#include <GH1D.h>

GH1D::GH1D() : TH1D() { 
}

GH1D::GH1D(const TH1D& rhs) : TH1D(rhs) { 
}

GH1D::GH1D(const char *name, int nbins, double low, double high) 
     : TH1D(name,name,nbins,low,high) {
}

GH1D::GH1D(const char *name,const char *title, int nbins, double low, double high) 
     : TH1D(name,title,nbins,low,high) {
}


GH1D::~GH1D() { }


void GH1D::Zoom(double low,double high) {
  if(low==low && high==high) {
    this->GetXaxis()->SetRangeUser(low,high);
  } else {
    this->GetXaxis()->UnZoom();
  }
}



int GH1D::Write(const char *name,Int_t option,Int_t bufsize) { 
  std::string hname = this->GetName();
  std::string temp_name = Form("__%s_temp__",this->GetName());
  this->SetName(temp_name.c_str());
  TH1D hist(hname.c_str(),this->GetTitle(),this->GetNbinsX(),
            this->GetXaxis()->GetBinLowEdge(1),this->GetXaxis()->GetBinUpEdge(this->GetNbinsX()));
  for(int i=0;i<=this->GetNbinsX()+1;i++) {
      hist.SetBinContent(i,this->GetBinContent(i));
  }
  hist.SetEntries(this->GetEntries()); 

  int result = hist.Write();
  hist.SetDirectory(0);
  hist.Delete();
  this->SetName(hname.c_str());
  return result;
}



