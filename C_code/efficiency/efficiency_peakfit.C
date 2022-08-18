
TH1 *NegRemove(TH1 *input){
  
  double nbins = input->GetNbinsX();
  double width = input->GetBinWidth(1);
  double lower = input->GetBinCenter(1);
  lower -= width/2;
  double upper = input->GetBinCenter(nbins);
  upper += width/2;
  TH1D *output = new TH1D(Form("%s_negremove",input->GetName()), Form("%s",input->GetTitle()), nbins,lower,upper);
  for(int i=1;i<=nbins;i++){
    if(input->GetBinContent(i)<0){
      //output->Fill(input->GetBinCenter(i),0);
      output->SetBinContent(i,0);
    }else{
      //output->Fill(input->GetBinCenter(i),input->GetBinContent(i));
      output->SetBinContent(i,input->GetBinContent(i));
    }
  }

  return output;
}






void SRM_PeakFit_xtal(int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp = (TH1D *)f->Get(Form("single_Cal_%i",xtal));

  int size = 9;
  double Energy[9] = {  105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower[9]  = {  102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double higher[9] = {  109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_xtal0_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("xtal%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    if(i<1) {
      gpeak0.push_back(GausFit(csp,lower[i],higher[i]));
    }else {
      gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    }
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    if(i<1){
      chi2 = ((GGaus *)gpeak0[i])->GetChisquare()/((GGaus *)gpeak0[i])->GetNDF();
      sum  = ((GGaus *)gpeak0[i])->GetSum();
    }else{
      chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
      sum  = ((GPeak *)gpeak0[i])->GetSum();
    }
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<Energy[i]<<"\t"<<((GGaus *)gpeak0[i])->GetSum()<<std::endl;
    else std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<((GGaus *)gpeak0[i])->GetSum()<<"\t";
    else std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}

void SRM_PeakFit_Sum(int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp;
  for(int i=4*xtal;i<4*xtal+4;i++){
    if(i==7 || i==9 || i==59) continue;
    TH1D *temp = (TH1D *)f->Get(Form("single_Cal_%i",i));
    if((i%4)==0){
      sp = (TH1D *)temp->Clone(Form("single_sum_%i",xtal));
    }else{
      sp->Add(temp,1);
    }
  }

  int size = 9;
  double Energy[9] = {  105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower[9]  = {  102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double higher[9] = {  109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_sum%i_peak%i_clone1",xtal,i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("Sum%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    if(i<1) {
      gpeak0.push_back(GausFit(csp,lower[i],higher[i]));
    }else {
      gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    }
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    if(i<1){
      chi2 = ((GGaus *)gpeak0[i])->GetChisquare()/((GGaus *)gpeak0[i])->GetNDF();
      sum  = ((GGaus *)gpeak0[i])->GetSum();
    }else{
      chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
      sum  = ((GPeak *)gpeak0[i])->GetSum();
    }
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<Energy[i]<<"\t"<<((GGaus *)gpeak0[i])->GetSum()<<std::endl;
    else std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<((GGaus *)gpeak0[i])->GetSum()<<"\t";
    else std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}

void SRM_PeakFit_Adb(int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp = (TH1D *)f->Get(Form("Addback_%i",xtal));

  int size = 9;
  double Energy[9] = {  105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower[9]  = {  102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double higher[9] = {  109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("Addback_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("Addback%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    if(i<1) {
      gpeak0.push_back(GausFit(csp,lower[i],higher[i]));
    }else {
      gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    }
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    if(i<1){
      chi2 = ((GGaus *)gpeak0[i])->GetChisquare()/((GGaus *)gpeak0[i])->GetNDF();
      sum  = ((GGaus *)gpeak0[i])->GetSum();
    }else{
      chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
      sum  = ((GPeak *)gpeak0[i])->GetSum();
    }
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<Energy[i]<<"\t"<<((GGaus *)gpeak0[i])->GetSum()<<std::endl;
    else std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<((GGaus *)gpeak0[i])->GetSum()<<"\t";
    else std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}



void SRM_PeakFit_AllSum(){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp;
  for(int i=0;i<64;i++){
    if(i==7 || i==9 || i==59) continue;
    TH1D *temp = (TH1D *)f->Get(Form("single_Cal_%i",i));
    if(i==0){
      sp = (TH1D *)temp->Clone("all_sum");
    }else{
      sp->Add(temp,1);
    }
  }

  int size = 9;
  double Energy[9] = {  105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower[9]  = {  102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double higher[9] = {  109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("all_sum_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("All Sum %.1fkeV bgsub negremove",Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    if(i<1) {
      gpeak0.push_back(GausFit(csp,lower[i],higher[i]));
    }else {
      gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    }
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    if(i<1){
      chi2 = ((GGaus *)gpeak0[i])->GetChisquare()/((GGaus *)gpeak0[i])->GetNDF();
      sum  = ((GGaus *)gpeak0[i])->GetSum();
    }else{
      chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
      sum  = ((GPeak *)gpeak0[i])->GetSum();
    }
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<Energy[i]<<"\t"<<((GGaus *)gpeak0[i])->GetSum()<<std::endl;
    else std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<((GGaus *)gpeak0[i])->GetSum()<<"\t";
    else std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}



void SRM_PeakFit_AllAdb(){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
  TH1D *sp;
  for(int i=0;i<16;i++){
    TH1D *temp = (TH1D *)f->Get(Form("Addback_%i",i));
    if(i==0){
      sp = (TH1D *)temp->Clone("all_addback");
    }else{
      sp->Add(temp,1);
    }
  }

  int size = 9;
  double Energy[9] = {  105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double lower[9]  = {  102.0, 118,   244,   587,   717,   865,   991,    1000,   1262};
  double higher[9] = {  109.5, 128,   252,   597,   730,   882,   1000.5, 1010,   1285};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("all_addback_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("All Addback %.1fkeV bgsub negremove",Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    if(i<1) {
      gpeak0.push_back(GausFit(csp,lower[i],higher[i]));
    }else {
      gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    }
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    if(i<1){
      chi2 = ((GGaus *)gpeak0[i])->GetChisquare()/((GGaus *)gpeak0[i])->GetNDF();
      sum  = ((GGaus *)gpeak0[i])->GetSum();
    }else{
      chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
      sum  = ((GPeak *)gpeak0[i])->GetSum();
    }
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<Energy[i]<<"\t"<<((GGaus *)gpeak0[i])->GetSum()<<std::endl;
    else std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    if(i<1) std::cout<<((GGaus *)gpeak0[i])->GetSum()<<"\t";
    else std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}
//============================= 56Co Peaks Fit =================================//
void Co_PeakFit_xtal(int xtal=0){
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TH1D *sp = (TH1D *)f->Get(Form("single_Cal_%i",xtal));

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  //double lower[12]  = {840,     1031,    1168,     1231,     1352,     1767,     2007,     2027,     2591,     3193,     3245,     3264};
  //double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};
  double lower[12]  = {840,     1031,    1166,     1229,     1352,     1760,     2007,     2027,     2590,     3193,     3245,     3264};
  double higher[12] = {855,     1046,    1183,     1247,     1369,     1779,     2022,     2042,     2606,     3210,     3262,     3282};


  TCanvas *c0 = new TCanvas;
  c0->Divide(3,4);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_xtal0_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("xtal%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}

void Co_PeakFit_Sum(int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TH1D *sp;
  for(int i=4*xtal;i<4*xtal+4;i++){
    if(i==7 || i==9 || i==59) continue;
    TH1D *temp = (TH1D *)f->Get(Form("single_Cal_%i",i));
    if((i%4)==0){
      sp = (TH1D *)temp->Clone(Form("single_sum_%i",xtal));
    }else{
      sp->Add(temp,1);
    }
  }

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  //double lower[12]  = {840,     1031,    1168,     1231,     1352,     1767,     2007,     2027,     2591,     3193,     3245,     3264};
  //double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};
  double lower[12]  = {840,     1031,    1166,     1229,     1352,     1760,     2007,     2027,     2590,     3193,     3245,     3264};
  double higher[12] = {855,     1046,    1183,     1247,     1369,     1779,     2022,     2042,     2606,     3210,     3262,     3282};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,4);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_sum%i_peak%i_clone1",xtal,i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("Sum%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}

void Co_PeakFit_Adb(int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TH1D *sp = (TH1D *)f->Get(Form("Addback_%i",xtal));

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  //double lower[12]  = {840,     1031,    1168,     1231,     1352,     1767,     2007,     2027,     2591,     3193,     3245,     3264};
  //double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};
  double lower[12]  = {840,     1031,    1166,     1229,     1352,     1760,     2007,     2027,     2590,     3193,     3245,     3264};
  double higher[12] = {855,     1046,    1183,     1247,     1369,     1779,     2022,     2042,     2606,     3210,     3262,     3282};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,4);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("Addback_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("Addback%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}


void Co_PeakFit_AllSum(){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TH1D *sp;
  for(int i=0;i<64;i++){
    if(i==7 || i==9 || i==59) continue;
    TH1D *temp = (TH1D *)f->Get(Form("single_Cal_%i",i));
    if(i==0){
      sp = (TH1D *)temp->Clone("all_sum");
    }else{
      sp->Add(temp,1);
    }
  }

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  //double lower[12]  = {840,     1031,    1168,     1231,     1352,     1767,     2007,     2027,     2591,     3193,     3245,     3264};
  //double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};
  double lower[12]  = {840,     1031,    1166,     1229,     1352,     1760,     2007,     2027,     2590,     3193,     3245,     3264};
  double higher[12] = {855,     1046,    1183,     1247,     1369,     1779,     2022,     2042,     2606,     3210,     3262,     3282};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,4);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("all_sum_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("All Sum %.1fkeV bgsub negremove",Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}


void Co_PeakFit_AllAdb(){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TH1D *sp;
  for(int i=0;i<16;i++){
    TH1D *temp = (TH1D *)f->Get(Form("Addback_%i",i));
    if(i==0){
      sp = (TH1D *)temp->Clone("all_addback");
    }else{
      sp->Add(temp,1);
    }
  }

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  //double lower[12]  = {840,     1031,    1168,     1231,     1352,     1767,     2007,     2027,     2591,     3193,     3245,     3264};
  //double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};
  double lower[12]  = {840,     1031,    1166,     1229,     1352,     1760,     2007,     2027,     2590,     3193,     3245,     3264};
  double higher[12] = {855,     1046,    1183,     1247,     1369,     1779,     2022,     2042,     2606,     3210,     3262,     3282};

  TCanvas *c0 = new TCanvas;
  c0->Divide(3,4);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  //new TCanvas;
  //sp->Draw();
  //bg->Draw("same");
  std::vector<TF1*> gpeak0;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("all_addback_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("All Addback %.1fkeV bgsub negremove",Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }

  std::cout<<"Energy"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"<<((GPeak *)gpeak0[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}









/*void Co_PeakFit_xtal(int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  TH1D *sp = (TH1D *)f->Get(Form("single_Cal_%i",xtal));

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  double lower[12]  = {840,     1031,    1168,     1231,     1352,     1767,     2007,     2027,     2591,     3193,     3245,     3264};
  double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};

  TCanvas *c0 = new TCanvas;
  TCanvas *c1 = new TCanvas;
  TCanvas *c2 = new TCanvas;
  c0->Divide(4,3);
  c1->Divide(4,3);
  c2->Divide(4,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  new TCanvas;
  sp->Draw();
  bg->Draw("same");
  std::vector<TF1*> gpeak0;
  std::vector<TF1*> gpeak1;
  std::vector<TF1*> gpeak2;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_xtal%i_peak%i_raw",xtal,i));
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("xtal%i %.1fkeV",xtal, Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }
  for(int i=0;i<size;i++){
    TH1D *csp = (TH1D *)sp->Clone(Form("single_xtal%i_peak%i_bgsub",xtal,i));
    csp->Add(bg,-1);
    csp->SetTitle(Form("xtal%i %.1fkeV bgsub",xtal, Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c1->cd(i+1);
    csp->Draw();
    gpeak1.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak1[i])->GetChisquare()/((GPeak *)gpeak1[i])->GetNDF();
    sum  = ((GPeak *)gpeak1[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c1->Update();
  }
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_xtal0_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("xtal%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c2->cd(i+1);
    csp->Draw();
    gpeak2.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak2[i])->GetChisquare()/((GPeak *)gpeak2[i])->GetNDF();
    sum  = ((GPeak *)gpeak2[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c2->Update();
  }

  std::cout<<"Energy"<<"\t"<<"sum"<<"\t"<<"bgsub"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"
             <<((GPeak *)gpeak0[i])->GetSum()<<"\t"
             <<((GPeak *)gpeak1[i])->GetSum()<<"\t"
             <<((GPeak *)gpeak2[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak1.size();i++){
    std::cout<<((GPeak *)gpeak1[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak2.size();i++){
    std::cout<<((GPeak *)gpeak2[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}




void Co_PeakFit(TH1D *sp=0, int xtal=0){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1143.root");
  if(sp == nullptr){
    bool flag = true;
    for(int i=xtal;i<xtal+4;i++){
      TH1D *temp = (TH1D *)f->Get(Form("single_Cal_%i",i));
      if(temp== nullptr) continue;
      if(flag){
        sp = (TH1D *)temp->Clone(Form("single_sum%i",xtal/4));
        sp->SetTitle(Form("Single Spectrum Sum Clover%i",xtal/4));
        flag = false;
      }else{
        sp->Add(temp,1);
      }
    }
  }

  int size = 12;
  double Energy[12] = {846.771, 1037.84, 1175.102, 1238.282, 1360.215, 1771.351, 2015.181, 2034.755, 2598.459, 3201.962, 3253.416, 3272.99};
  //double lower[12]  = {840,     1031,    1167,     1232,     1352,     1766,     2007,     2027,     2591,     3193,     3245,     3265};
  //double higher[12] = {855,     1043,    1182,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};
  double lower[12]  = {840,     1031,    1168,     1230.5,     1352,     1766,     2007,     2027,     2591,     3193,     3244,     3265};
  double higher[12] = {855,     1043,    1181,     1243,     1365,     1778,     2022,     2041,     2606,     3210,     3262,     3282};

  TCanvas *c0 = new TCanvas;
  TCanvas *c1 = new TCanvas;
  TCanvas *c2 = new TCanvas;
  c0->Divide(4,3);
  c1->Divide(4,3);
  c2->Divide(4,3);
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  new TCanvas;
  sp->Draw();
  bg->Draw("same");
  std::vector<TF1*> gpeak0;
  std::vector<TF1*> gpeak1;
  std::vector<TF1*> gpeak2;
  TText *txt[3];
  double chi2, sum;
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_xtal%i_peak%i_raw",xtal,i));
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("xtal%i %.1fkeV",xtal, Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c0->cd(i+1);
    csp->Draw();
    gpeak0.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak0[i])->GetChisquare()/((GPeak *)gpeak0[i])->GetNDF();
    sum  = ((GPeak *)gpeak0[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c0->Update();
  }
  for(int i=0;i<size;i++){
    TH1D *csp = (TH1D *)sp->Clone(Form("single_xtal%i_peak%i_bgsub",xtal,i));
    csp->Add(bg,-1);
    csp->SetTitle(Form("xtal%i %.1fkeV bgsub",xtal, Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c1->cd(i+1);
    csp->Draw();
    gpeak1.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak1[i])->GetChisquare()/((GPeak *)gpeak1[i])->GetNDF();
    sum  = ((GPeak *)gpeak1[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c1->Update();
  }
  for(int i=0;i<size;i++){
    TH1D *clsp = (TH1D *)sp->Clone(Form("single_xtal0_peak%i_clone1",i));
    clsp->Add(bg,-1);
    TH1D *csp = (TH1D *)NegRemove(clsp);
    csp->SetTitle(Form("xtal%i %.1fkeV bgsub negremove",xtal,Energy[i]));
    csp->GetXaxis()->SetRangeUser(lower[i]-20, higher[i]+20);
    csp->GetListOfFunctions()->Clear();
    c2->cd(i+1);
    csp->Draw();
    gpeak2.push_back(PhotoPeakFit(csp,lower[i],higher[i])); 
    for(int txtc=0;txtc<3;txtc++){
      txt[txtc] = new TText();
      txt[txtc] -> SetNDC();
      txt[txtc] -> SetTextFont(1);
      txt[txtc] -> SetTextColor(1);
      txt[txtc] -> SetTextSize(0.05);
    }
    double intg = csp->Integral(csp->GetXaxis()->FindBin(lower[i]),
                                csp->GetXaxis()->FindBin(higher[i]));
    chi2 = ((GPeak *)gpeak2[i])->GetChisquare()/((GPeak *)gpeak2[i])->GetNDF();
    sum  = ((GPeak *)gpeak2[i])->GetSum();
    txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
    txt[1] -> DrawText(0.1,0.6, Form("integral = %.3f",intg));
    txt[2] -> DrawText(0.1,0.5, Form("sum  = %.3f",sum));
    c2->Update();
  }

  std::cout<<"Energy"<<"\t"<<"sum"<<"\t"<<"bgsub"<<"\t"<<"NegRemove"<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t"
             <<((GPeak *)gpeak0[i])->GetSum()<<"\t"
             <<((GPeak *)gpeak1[i])->GetSum()<<"\t"
             <<((GPeak *)gpeak2[i])->GetSum()<<std::endl;
  }
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  for(int i=0;i<gpeak0.size();i++){
    std::cout<<Energy[i]<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak0.size();i++){
    std::cout<<((GPeak *)gpeak0[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak1.size();i++){
    std::cout<<((GPeak *)gpeak1[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  for(int i=0;i<gpeak2.size();i++){
    std::cout<<((GPeak *)gpeak2[i])->GetSum()<<"\t";
  }
  std::cout<<std::endl;
  return;
}*/


//============================ Single Peak Fit =======================//
void SRMSinglePeakFit(TH1D *sp=0, int en=2, Option_t *opt="sum"){

  double energy[10] = {86.5, 105.3, 123.1, 247.7, 591.8, 723.3, 873.2, 996.3,  1004.7, 1274.5};
  double emit[4]    = {868008.225421772, 604533.07196845, 27253844.9047775, 4618975.33299466};
  double lower[4] = {82, 103.5, 118, 244};
  double upper[4] = {90, 108,   128, 252};

  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");

  if(sp==nullptr){
    TFile *f = TFile::Open("/home/zhu/packages/BCSSort/spectrum/list_output1142.root");
    TH1D *Addback_0 = (TH1D *)f->Get("Addback_0");
    TH1D *single_Cal_0 = (TH1D *)f->Get("single_Cal_0");
    TH1D *sum0; 
    for(int i=0;i<4;i++){
      TH1D *temp = (TH1D *)f->Get(Form("single_Cal_%i",i));
      if(i==0){
        sum0 = (TH1D *)temp->Clone("sum0_clone");
      }else{
        sum0->Add(temp,1);
      }
    }
    if(sopt.Contains("sum")){
      sp = (TH1D *)sum0->Clone("sim0");
      sp->SetTitle("sum0");
    } 
    if(sopt.Contains("addb")){
      sp = (TH1D *)Addback_0->Clone("addback0");
      sp->SetTitle("addback0");
    } 
  }
  
  TSpectrum s;
  TH1D *bg = (TH1D *)s.Background(sp,10);
  TH1D *csp = (TH1D *)sp->Clone(Form("%s_bgsub",sp->GetName()));
  csp->Add(bg,-1);
  TH1D *ncsp = (TH1D *)NegRemove(csp);
  ncsp->GetXaxis()->SetRangeUser(lower[en]-20,upper[en]+20);
  new TCanvas;
  ncsp->Draw();
  gPad->SetLogy();  

  std::vector<GPeak *> gpeak;
  std::vector<GPeak *> gpeak1;
  TCanvas *c = new TCanvas;
  c->Divide(4,4);
  TCanvas *c1 = new TCanvas;
  c1->Divide(4,4);
  TText *txt[2];
  for(int i=0;i<2;i++){
    txt[i] = new TText;
    txt[i] -> SetNDC();
    txt[i] -> SetTextFont(1);
    txt[i] -> SetTextColor(1);
    txt[i] -> SetTextSize(0.05);
  }
  TText *txt1[2];
  for(int i=0;i<2;i++){
    txt1[i] = new TText;
    txt1[i] -> SetNDC();
    txt1[i] -> SetTextFont(1);
    txt1[i] -> SetTextColor(1);
    txt1[i] -> SetTextSize(0.05);
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      double xlower = lower[en]-0.5*j;
      double xupper = upper[en]+0.5*i;
      TH1D *csp = (TH1D *)sp->Clone(Form("%s_%ibgsub",sp->GetName(),4*i+j));
      csp->Add(bg,-1);
      TH1D *ncsp = (TH1D *)NegRemove(csp);
      ncsp->GetXaxis()->SetRangeUser(lower[en]-20,upper[en]+20);
      ncsp->GetListOfFunctions()->Clear();
      ncsp->SetTitle(Form("%s %.1fkeV [%.1f,%.1f]", sp->GetTitle(),energy[en], xlower, xupper));
      TH1D *cncsp = (TH1D *)ncsp->Clone(Form("%s_clone",ncsp->GetName()));
      cncsp->GetXaxis()->SetRangeUser(lower[en]-20,upper[en]+20);
      cncsp->GetListOfFunctions()->Clear();
      c->cd(1+i*4+j);
      ncsp->Draw();
      gPad->SetLogy();
      gpeak.push_back(PhotoPeakFit(ncsp,xlower,xupper));
      double chi2 = ((GPeak *)gpeak[i+j])->GetChisquare()/((GPeak *)gpeak[i+j])->GetNDF();
      double sum  = ((GPeak *)gpeak[i+j])->GetSum();
      txt[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
      txt[1] -> DrawText(0.1,0.6, Form("sum  = %.3f",sum));
      c1->cd(1+i*4+j);
      cncsp->GetYaxis()->SetRangeUser(10,1000);
      cncsp->Draw();
      gPad->SetLogy();
      gpeak1.push_back(PhotoPeakFit(cncsp,xlower,xupper));
      chi2 = ((GPeak *)gpeak[i+j])->GetChisquare()/((GPeak *)gpeak[i+j])->GetNDF();
      sum  = ((GPeak *)gpeak[i+j])->GetSum();
      txt1[0] -> DrawText(0.1,0.7, Form("chi2 = %.3f",chi2));
      txt1[1] -> DrawText(0.1,0.6, Form("sum  = %.3f",sum));
      c->Update(); 
      c1->Update(); 
    }
  }

  std::cout<<"lower" << "\t" << "upper" << "\t" << "sum" << "\t" <<"eff" << "\t"<<"effdif"<<std::endl;
  double eff0 = gpeak[0]->GetSum()/emit[en];
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      double xlower = lower[en]-0.5*j;
      double xupper = upper[en]+0.5*i;
      double eff = gpeak[i+j]->GetSum()/emit[en];
      double dif = eff-eff0;
      std::cout<<xlower<<"\t"<<xupper<<"\t"<<gpeak[i+j]->GetSum()<<"\t"<<eff<<"\t"<<dif<<std::endl;
    }
  } 

  return;




}


