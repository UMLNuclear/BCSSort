


// Run in GRUTinizer (because of "GGasu")//

void tofcorrection(){

  TFile *f;
  double sigma = 20;  
  ofstream outfile;
  outfile.open("ctof_peaks.dat");

  for(int x=1004;x<1140;x++){
    f = TFile::Open(Form("root_file/tof_correction/ctof/ctof%i.root",x));
    if(f==nullptr){
      printf("ctof%i.root doesn't exist\n",x);
      continue;
    }
    
    TH2D *ctof = (TH2D *)f->Get(Form("ctof%i",x));
    if((ctof->Integral())<10){
      printf("ctof%i.root is empty\n",x);
      continue;
    }

    TH1D *pctof = ctof->ProjectionY(Form("ctof_py_%i",x));
    TSpectrum s;
    s.Search(pctof,sigma,"",0.25);
    vector<double> cen;
    for(int i=0;i<s.GetNPeaks();i++){
      GGaus *peak = GausFit(pctof, s.GetPositionX()[i]-2*sigma, s.GetPositionX()[i]+2*sigma);
      //pctof->GetXaxis()->SetRangeUser(s.GetPositionX()[i]-2*sigma, s.GetPositionX()[i]+2*sigma);
      //pctof->Fit("gaus");
      //TF1 *g = (TF1 *)pctof->GetListOfFunctions()->FindObject("gaus");
      //cen.push_back(g->GetParameter(1));
      cen.push_back(peak->GetCentroid());
      
    }
    sort(cen.begin(),cen.end());
    outfile<< setw(5)<<x;
    for(int i=0;i<cen.size();i++){
      outfile<<"\t"<<std::setw(10) <<cen[i] << "\t";      
    }
    outfile<<endl;

  }

  outfile.close();
}
