

void PIDAnaly(){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/PID_gateGam_exdTOF_corI2.root");
  TFile *cutf = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut.root");
  TCutG *na = (TCutG *)cutf->Get("Na");
  na->SetName("na");
  TCutG *ne = (TCutG *)cutf->Get("Ne");
  ne->SetName("ne");
  TCutG *nacor = new TCutG("nacor",5);
  nacor->SetName("nacor");
  nacor->SetPoint(0,15200,6300);
  nacor->SetPoint(1,15200,7300);
  nacor->SetPoint(2,19500,6600);
  nacor->SetPoint(3,19500,5600);
  nacor->SetPoint(4,15200,6300);
  nacor->SetLineColor(kBlack);
  nacor->SetLineWidth(2);
  TCutG *necor = new TCutG("necor",5);
  necor->SetName("necor");
  necor->SetPoint(0,13500,5400);
  necor->SetPoint(1,13500,6400);
  necor->SetPoint(2,19100,5500);
  necor->SetPoint(3,19100,4500);
  necor->SetPoint(4,13500,5400);
  necor->SetLineColor(kBlack);
  necor->SetLineWidth(2);
  TCanvas *c1 = new TCanvas;
  c1->Divide(1,3);
  TCanvas *c2 = new TCanvas;
  c2->Divide(1,3);
  TCanvas *c3 = new TCanvas;
  c3->Divide(1,3);
  TCanvas *c4 = new TCanvas;
  c4->Divide(1,3);
  TCanvas *c5 = new TCanvas;
  c5->Divide(1,3);
  TCanvas *c6 = new TCanvas;
  c6->Divide(1,3);
  TCanvas *c7 = new TCanvas;
  c7->Divide(1,3);

  //=============== Gate 150keV ===========//
  TH2D *pid = (TH2D *)f->Get("pid_100ms_885W");
  pid->SetTitle("PID gated[882,888]keV within 100ms");
  TH2D *pid_BGR = (TH2D *)f->Get("pid_100ms_BGR885W");
  pid_BGR->SetTitle("PID gated[900,906]keV within 100ms");
  TH2D *dtsingles = (TH2D *)f->Get("dt_singles");
  TH1D *singles = (TH1D *)dtsingles->ProjectionY("singles_na",1,100);
  singles->GetXaxis()->SetRangeUser(880,910);
  singles->SetTitle("singles within 100ms");
  c1->cd(1);
  pid->Draw("colz");
  na->Draw("same");
  ne->Draw("same");
  c1->cd(2);
  pid_BGR->Draw("colz");
  na->Draw("same");
  ne->Draw("same");
  c1->cd(3);
  singles->Draw();
  TLine *l1 = new TLine(882,0,882,singles->GetMaximum());
  l1->SetLineColor(kRed); 
  l1->SetLineWidth(2);
  l1->Draw("same");   
  TLine *l2 = new TLine(888,0,888,singles->GetMaximum());
  l2->SetLineColor(kRed); 
  l2->SetLineWidth(2);
  l2->Draw("same");   
  TLine *l3 = new TLine(900,0,900,singles->GetMaximum());
  l3->SetLineColor(kBlue); 
  l3->SetLineWidth(2);
  l3->Draw("same");   
  TLine *l4 = new TLine(906,0,906,singles->GetMaximum());
  l4->SetLineColor(kBlue); 
  l4->SetLineWidth(2);
  l4->Draw("same");   

  TH1D *tofna = pid->ProjectionX("tof_na", 1,300,"[na]");  
  TH1D *tofbgna = pid_BGR->ProjectionX("tof_bg_na", 1,300,"[na]");  
  TH1D *ctofna = (TH1D *)tofna->Clone(Form("%s_clone",tofna->GetName()));
  ctofna->Add(tofbgna,-1);
  c2->cd(1);
  tofna->Draw();
  c2->cd(2);
  tofbgna->Draw();
  c2->cd(3);
  tofna->Draw();
  ctofna->SetLineColor(kRed);
  ctofna->Draw("same");

  TH1D *tofne = pid->ProjectionX("tof_ne", 1,300,"[ne]");  
  TH1D *tofbgne = pid_BGR->ProjectionX("tof_bg_ne", 1,300,"[ne]");  
  TH1D *ctofne = (TH1D *)tofne->Clone(Form("%s_clone",tofne->GetName()));
  ctofne->Add(tofbgne,-1);
  c3->cd(1);
  tofne->Draw();
  c3->cd(2);
  tofbgne->Draw();
  c3->cd(3);
  ctofne->Draw();

  //=============== Gate 150keV ===========//
  TH2D *pidne = (TH2D *)f->Get("pid_50ms_150W");
  pidne->SetTitle("PID gated[146,152]keV within 50ms");
  TH2D *pidne_BGR = (TH2D *)f->Get("pid_50ms_BGR150W");
  pidne_BGR->SetTitle("PID gated[156,162]keV within 50ms");
  TH1D *singlesne = (TH1D *)dtsingles->ProjectionY("singles_ne",1,50);
  singlesne->GetXaxis()->SetRangeUser(140,165);
  singlesne->SetTitle("singles within 50ms");
  c4->cd(1);
  pidne->Draw("colz");
  na->Draw("same");
  ne->Draw("same");
  c4->cd(2);
  pidne_BGR->Draw("colz");
  na->Draw("same");
  ne->Draw("same");
  c4->cd(3);
  singlesne->Draw();
  TLine *l5 = new TLine(146,0,146,singlesne->GetMaximum());
  l5->SetLineColor(kRed); 
  l5->SetLineWidth(2);
  l5->Draw("same");   
  TLine *l6 = new TLine(152,0,152,singlesne->GetMaximum());
  l6->SetLineColor(kRed); 
  l6->SetLineWidth(2);
  l6->Draw("same");   
  TLine *l7 = new TLine(156,0,156,singlesne->GetMaximum());
  l7->SetLineColor(kBlue); 
  l7->SetLineWidth(2);
  l7->Draw("same");   
  TLine *l8 = new TLine(162,0,162,singlesne->GetMaximum());
  l8->SetLineColor(kBlue); 
  l8->SetLineWidth(2);
  l8->Draw("same");   

  TH1D *tofNEna = pidne->ProjectionX("tofNE_na", 1,300,"[na]");  
  TH1D *tofNEbgna = pidne_BGR->ProjectionX("tofNE_bg_na", 1,300,"[na]");  
  TH1D *ctofNEna = (TH1D *)tofNEna->Clone(Form("%s_clone",tofNEna->GetName()));
  ctofNEna->Add(tofNEbgna,-1);
  c5->cd(1);
  tofNEna->Draw();
  c5->cd(2);
  tofNEbgna->Draw();
  c5->cd(3);
  ctofNEna->Draw();

  TH1D *tofNEne = pidne->ProjectionX("tofNE_ne", 1,300,"[ne]");  
  TH1D *tofNEbgne = pidne_BGR->ProjectionX("tofNE_bg_ne", 1,300,"[ne]");  
  TH1D *ctofNEne = (TH1D *)tofNEne->Clone(Form("%s_clone",tofNEne->GetName()));
  ctofNEne->Add(tofNEbgne,-1);
  c6->cd(1);
  tofNEne->Draw();
  c6->cd(2);
  tofNEbgne->Draw();
  c6->cd(3);
  tofNEne->Draw();
  ctofNEne->SetLineColor(kRed);
  ctofNEne->Draw("same");

  //============ Compare pid gate 150keV before and after TOF correction(gate I2_pos or not)

  TFile *f1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/PID_gate.root");
  TH2D *pidne_bef = (TH2D *)f1->Get("pid_50ms_150");
  pidne_bef->SetTitle("PID_beforeTOFexf gated[146,152]keV within 50ms");
  TH2D *pidne_bef_BGR = (TH2D *)f1->Get("pid_50ms_BGR150");
  pidne_bef_BGR->SetTitle("PID_beforeTOFexd gated[156,162]keV within 50ms");
  TH1D *tofbefNEne = pidne_bef->ProjectionX("tofbefNE_ne", 1,300,"[ne]");  
  TH1D *tofbefNEbgne = pidne_bef_BGR->ProjectionX("tofbefNE_bg_ne", 1,300,"[ne]");  
  TH1D *ctofbefNEne = (TH1D *)tofbefNEne->Clone(Form("%s_clone",tofbefNEne->GetName()));
  ctofbefNEne->Add(tofbefNEbgne,-1);

  TH2D *pidcorne = (TH2D *)f->Get("pid_50ms_150_cor");
  pidcorne->SetTitle("PIDcor gated[146,152]keV within 50ms");
  TH2D *pidcorne_BGR = (TH2D *)f->Get("pid_50ms_BGR150_cor");
  pidcorne_BGR->SetTitle("PIDcor gated[156,162]keV within 50ms");
  TH1D *tofcorNEne = pidcorne->ProjectionX("tofcorNE_ne", 1,300,"[necor]");  
  TH1D *tofcorNEbgne = pidcorne_BGR->ProjectionX("tofcorNE_bg_ne", 1,300,"[necor]");  
  TH1D *ctofcorNEne = (TH1D *)tofcorNEne->Clone(Form("%s_clone",tofcorNEne->GetName()));
  ctofcorNEne->Add(tofcorNEbgne,-1);

  c7->cd(1);
  tofbefNEne->GetXaxis()->SetRangeUser(10000,16000);
  tofbefNEne->Draw();
  ctofbefNEne->SetLineColor(kRed);
  ctofbefNEne->Draw("same");
  c7->cd(2);
  tofNEne->GetXaxis()->SetRangeUser(10000,16000);
  tofNEne->Draw();
  ctofNEne->Draw("same");
  c7->cd(3); 
  tofcorNEne->Draw();
  tofcorNEne->GetXaxis()->SetRangeUser(13000,19000);
  ctofcorNEne->SetLineColor(kRed);
  ctofcorNEne->Draw("same");


  return;

}



void NasingleDraw(){

  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdCor_slide523.root");
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  double eadd = 882;
  double esub = 900;
  double gsize = 6;
  double time = 100.;
  double tbin = time*1000;
  int rbin = 2000;
  TH2D *pidna = (TH2D *)f->Get("PID_Na");

  TCanvas *c1 = new TCanvas;
  TCanvas *c2 = new TCanvas;
  TCanvas *c3 = new TCanvas;
  TCanvas *c4 = new TCanvas;
  TCanvas *c5 = new TCanvas;
  TCanvas *c6 = new TCanvas;
  c1->Divide(5,2);
  c2->Divide(5,2);
  c3->Divide(1,10);
  c4->Divide(5,2);
  c5->Divide(5,2);
  c6->Divide(5,2);

  TGraphErrors *gr = new TGraphErrors(8);
  gr->SetTitle(Form("885keV, first %.1fms",time));
  TGraphErrors *gr1 = new TGraphErrors(12);
  gr1->SetTitle(Form("Chi2 885keV, first %.1fms",time));
  for(int i=0;i<10;i++){
    TH2D *pidcut = (TH2D *)f->Get(Form("PID_recna%i",i));
    pidcut->SetTitle(Form("PID_Na_cut%i",i));
    c1->cd(i+1);
    pidna->Draw("colz");
    pidcut->Draw("same");

    TH2D *i2 = (TH2D *)f->Get(Form("I2_TOF_recna%i",i));
    i2->SetTitle(Form("I2_Na_cut%i",i));
    i2->GetYaxis()->SetRangeUser(7000,17000);
    //i2->GetXaxis()->SetRangeUser(12500+i*230,13200+i*230);
    c2->cd(i+1);
    i2->Draw("colz");

    TH2D *dtsingles = (TH2D *)f->Get(Form("dt_singles_recna%i",i));
    TH1D *singles = dtsingles->ProjectionY(Form("singles_100ms_%i",i),1,tbin);
    c3->cd(i+1);
    singles->Draw();

    TH1D *csingles = (TH1D *)singles->Clone(Form("%s_clone",singles->GetName())); 
    csingles->GetXaxis()->SetRangeUser(eadd-gsize,esub+2*gsize);
    c5->cd(i+1);
    csingles->Draw();
    TLine *l1 = new TLine(eadd,0,eadd,singles->GetMaximum());
    l1->SetLineColor(kRed); 
    l1->SetLineWidth(2);
    l1->Draw("same");   
    TLine *l2 = new TLine(eadd+gsize,0,eadd+gsize,singles->GetMaximum());
    l2->SetLineColor(kRed); 
    l2->SetLineWidth(2);
    l2->Draw("same");   
    TLine *l3 = new TLine(esub,0,esub,singles->GetMaximum());
    l3->SetLineColor(kBlue); 
    l3->SetLineWidth(2);
    l3->Draw("same");   
    TLine *l4 = new TLine(esub+gsize,0,esub+gsize,singles->GetMaximum());
    l4->SetLineColor(kBlue); 
    l4->SetLineWidth(2);
    l4->Draw("same");   

    TH1D *gadd = dtsingles->ProjectionX(Form("dt_885keV_%i",i),eadd,eadd+gsize);
    TH1D *gsub = dtsingles->ProjectionX(Form("dt_900keV_%i",i),esub,esub+gsize);

    //TH1D *cgadd = (TH1D *)gadd->Clone(Form("%s_clone",gadd->GetName()));
    //cgadd->Add(gsub,-1);
    //cgadd->Rebin(2000);
    //cgadd->GetXaxis()->SetRangeUser(time-20,time);
    //c6->cd(i+1);
    //cgadd->Draw();
    //TF1 *fbg = new TF1("fbg","[0]");
    //cgadd->Fit(fbg);

    gadd->Add(gsub,-1);
    gadd->Rebin(rbin);
    gadd->GetXaxis()->SetRangeUser(0,time);
    c4->cd(i+1);
    gadd->Draw();
    TF1 *fun = new TF1("fun","[0]*exp(-0.693/[1]*x)+[2]",0,100);
    fun->SetParameters(gadd->GetMaximum(),10,0.01);
    //fun->FixParameter(2,fbg->GetParameter(0));
    gadd->Fit(fun);
    TText *chi2 = new TText();
    chi2 -> SetNDC();
    chi2 -> SetTextFont(1);
    chi2 -> SetTextColor(1);
    chi2 -> SetTextSize(0.05);
    chi2 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl = new TText();
    hl -> SetNDC();
    hl -> SetTextFont(1);
    hl -> SetTextColor(1);
    hl -> SetTextSize(0.05);
    hl -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1),fun->GetParError(1)));
    gr->SetPoint(i,i+1,fun->GetParameter(1));
    gr->SetPointError(i,0,fun->GetParError(1)); 
    gr1->SetPoint(i,i+1,fun->GetChisquare()/fun->GetNDF());


  }
  new TCanvas;
  gr->Draw("A*");
  new TCanvas;
  gr1->Draw("A*");


}




void CutDraw(){

  TCutG *cutg;
  double lna = 700;
  double wna = 800;
  double horna = 230;
  double verna = 40;
  for(int i=0;i<10;i++){
    double start[2] = {12500,6200};
    start[0] += i*horna;
    start[1] -= i*verna;
    cutg = new TCutG(Form("recna%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]+wna);
    cutg->SetPoint(3,start[0],start[1]+wna);
    cutg->SetPoint(4,start[0],start[1]);
    cutg->SetLineColor(i+1);
    cutg->SetLineWidth(2);
    cutg->Draw("same");
  }

  TCutG *NaCut = new TCutG("NaCut",5);
  NaCut->SetPoint(0,12500,6150);
  NaCut->SetPoint(1,15500,5700);
  NaCut->SetPoint(2,15500,6700);
  NaCut->SetPoint(3,12500,7200);
  NaCut->SetPoint(4,12500,6150);
  NaCut->SetLineWidth(4);
  NaCut->Draw("same");

}





void NesingleDraw(){

  // bin count;
  double eadd = 146;
  //double esub = 156;
  double esub = 152;
  double gsize = 6;
  double time = 50; //unit:ms
  double starttime = 0; // unit:ms
  double startbin = starttime*100; // unit:us
  double tbin = time*100;  // unit:us
  int rbin = 50;


  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdUnCor_slide547_dtbin1us.root");
  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  TCanvas *c1 = new TCanvas;
  TCanvas *c2 = new TCanvas;
  TCanvas *c3 = new TCanvas;
  TCanvas *c4 = new TCanvas;
  TCanvas *c5 = new TCanvas;
  TCanvas *c6 = new TCanvas;
  //TCanvas *c7 = new TCanvas;
  c1->Divide(5,4);
  c2->Divide(5,4);
  c3->Divide(2,10);
  c4->Divide(5,4);
  c5->Divide(5,4);
  c6->Divide(5,4);
  //c7->Divide(4,4);
  TGraphErrors *gr = new TGraphErrors(20);
  gr->SetTitle(Form("150keV, first %.1fms",time));
  TGraphErrors *gr1 = new TGraphErrors(20);
  gr1->SetTitle(Form("Chi2 150keV, first %.1fms",time));

  for(int i=0;i<20;i++){

    TH2D *PID = (TH2D *)f->Get("pid_before");
    PID->GetXaxis()->SetRangeUser(9000,16000);
    PID->GetYaxis()->SetRangeUser(4700,6500);
    PID->SetTitle(Form("PID_Ne_cut%i",i));
    c1->cd(i+1);
    PID->Draw("colz");
    TH2D *pid = (TH2D *)f->Get(Form("PID_rec%i",i));
    TH2D *cpid = (TH2D *)pid->Clone(Form("%s_clone",pid->GetName()));
    cpid->Draw("same"); 

    TH2D *i2p = (TH2D *)f->Get(Form("I2_TOF_rec%i",i));
    i2p->SetTitle(Form("I2_TOF_Ne_cut%i",i));
    c2->cd(i+1);
    i2p->Draw("colz");


    TH2D *dtsingles = (TH2D *)f->Get(Form("dt_singles_rec%i",i));
    TH1D *singles = dtsingles->ProjectionY(Form("singles_rec%i",i),startbin,tbin);
    TH1D *gadd = dtsingles->ProjectionX(Form("dt150_rec%i",i),eadd,eadd+gsize);
    TH1D *gsub = dtsingles->ProjectionX(Form("dtbg_rec%i",i),esub, esub+36);
    singles->SetStats(0);
    c3->cd(i+1);
    singles->Draw();
    TH1D *cgsub = (TH1D *)gsub->Clone(Form("%s_clone",gsub->GetName()));
    cgsub->Rebin(rbin);
    cgsub->GetXaxis()->SetRangeUser(starttime,time);
    c6->cd(i+1);
    cgsub->Draw();

    TH1D *cgadd = (TH1D *)gadd->Clone(Form("%s_clone",gadd->GetName()));
    cgadd->Add(gsub,-1/6.);
    cgadd->Rebin(rbin);
    cgadd->GetXaxis()->SetRangeUser(time,100);
    TF1 *fbg = new TF1("fbg","[0]");
    cgadd->Fit(fbg);
    
    TH1D *csingles = (TH1D *)singles->Clone(Form("%s_clone",singles->GetName()));  
    csingles->GetXaxis()->SetRangeUser(eadd-gsize,esub+2.*gsize);
    c4->cd(i+1);
    csingles->Draw();
    TLine *l1 = new TLine(eadd,0,eadd,csingles->GetMaximum());
    l1->SetLineColor(kRed); 
    l1->SetLineWidth(2);
    l1->Draw("same");   
    TLine *l2 = new TLine(eadd+gsize,0,eadd+gsize,csingles->GetMaximum());
    l2->SetLineColor(kRed); 
    l2->SetLineWidth(2);
    l2->Draw("same");   
    TLine *l3 = new TLine(esub,0,esub,csingles->GetMaximum());
    l3->SetLineColor(kBlue); 
    l3->SetLineWidth(2);
    l3->Draw("same");   
    TLine *l4 = new TLine(esub+gsize,0,esub+gsize*6,csingles->GetMaximum());
    l4->SetLineColor(kBlue); 
    l4->SetLineWidth(2);
    l4->Draw("same");   

    gadd->Add(gsub,-1/6.);
    gadd->Rebin(rbin);
    gadd->GetXaxis()->SetRangeUser(starttime,time);
    c5->cd(i+1);
    gadd->Draw();
    TF1 *fun = new TF1(Form("fun_%i",i),"[0]*exp(-0.693/[1]*x)+[2]",0,50);
    fun->SetParLimits(2,0,20);
    fun->SetParameters(gadd->GetMaximum(),5,0.01);
    //fun->FixParameter(2,fbg->GetParameter(0));
    printf("\n\ncut_%i",i);
    gadd->Fit(fun);
    TText *chi2 = new TText();
    chi2 -> SetNDC();
    chi2 -> SetTextFont(1);
    chi2 -> SetTextColor(1);
    chi2 -> SetTextSize(0.05);
    chi2 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl = new TText();
    hl -> SetNDC();
    hl -> SetTextFont(1);
    hl -> SetTextColor(1);
    hl -> SetTextSize(0.05);
    hl -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));
    gr->SetPoint(i,i+1,fun->GetParameter(1));
    gr->SetPointError(i,0,fun->GetParError(1));
    gr1->SetPoint(i,i+1,fun->GetChisquare()/fun->GetNDF());


  }
  new TCanvas;
  gr->Draw("A*");
  new TCanvas;
  gr1->Draw("A*");

}



void SearchSubGate(){

  double eadd = 146;
  double esub = 153;
  double gsize = 6;
  double time = 50;

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdUnCor_slide535.root");
  TCanvas *c;
  TCanvas *c100 = new TCanvas;
  c100->Divide(5,2);
  TGraphErrors *gr;
  for(int j=0;j<10;j++){
    gr = new TGraphErrors(8);
    esub = 153 + j;
    c = new TCanvas;
    c->Divide(4,2);
    for(int i=0;i<8;i++){
      TH2D *dtsingles = (TH2D *)f->Get(Form("dt_singles_rec%i",i));
      TH1D *gsub = dtsingles->ProjectionX(Form("dtbg_%.1fkeV_rec%i",esub,i),esub, esub+gsize);
      gsub->SetTitle(Form("dtbg @ %.2f rec%i",esub,i));
      gsub->Rebin(2);
      gsub->GetXaxis()->SetRangeUser(1,time*2);
      c->cd(i+1);
      gsub->Draw();
      TF1 *fun = new TF1("fun","[0]*exp(-0.693/[1]*x)+[2]");
      fun->SetParameters(gsub->GetMaximum(),50,gsub->GetMinimum());
      gsub->Fit(fun);
      gr->SetPoint(i,i+1,fun->GetParameter(1));
      gr->SetPointError(i,0,fun->GetParError(1));
    }
    c100->cd(j+1);
    gr->Draw("A*");
  }
  
  c = new TCanvas;
  c->Divide(4,2);
  gr = new TGraphErrors(8);
  for(int i=0;i<8;i++){
    TH2D *dtsingles = (TH2D *)f->Get(Form("dt_singles_rec%i",i));
    TH1D *gadd = dtsingles->ProjectionX(Form("dtadd_%.1fkeV_rec%i",eadd,i),eadd, eadd+gsize);
    gadd->Rebin(2);
    gadd->GetXaxis()->SetRangeUser(1,time*2);
    c->cd(i+1);
    gadd->Draw();
    TF1 *fun = new TF1("fun","[0]*exp(-0.693/[1]*x)+[2]");
    fun->SetParameters(gadd->GetMaximum(),6,gadd->GetMinimum());
    gadd->Fit(fun);
    gr->SetPoint(i,i+1,fun->GetParameter(1));
    gr->SetPointError(i,0,fun->GetParError(1));
  }
  new TCanvas;
  gr->Draw("A*");

  return;
}


void NeDifBinWidth() {

  // bin count;
  double eadd = 146;
  double esub = 154;
  double gsize = 6;
  double time = 50; //unit:ms
  double starttime = 0.01; // unit:ms
  double startbin = 1*1000; // unit:us
  double tbin = time*1000;  // unit:us
  int rbin = 10;

  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdUnCor_slide535.root");
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  
  TCanvas *c1 = new TCanvas;
  c1->Divide(4,5);
  TCanvas *c2 = new TCanvas;
  c2->Divide(4,5);

  TGraphErrors *gr = new TGraphErrors(20);
  gr->SetTitle(Form("150keV, first %.1fms",time));
  TGraphErrors *gr1 = new TGraphErrors(20);
  gr1->SetTitle(Form("Chi2 150keV, first %.1fms",time));

  TGraphErrors *gr2 = new TGraphErrors(20);
  gr2->SetTitle(Form("150keV, first %.1fms(no 1st bin)",time));
  TGraphErrors *gr21 = new TGraphErrors(20);
  gr21->SetTitle(Form("Chi2 150keV, first %.1fms(no 1st bin)",time));
  
  int i = 5;
  TH2D *dtsingles = (TH2D *)f->Get(Form("dt_singles_rec%i",i));
  TH1D *dt   = dtsingles->ProjectionX(Form("dt_rec%i",i),eadd, eadd+gsize); 
  TH1D *gsub = dtsingles->ProjectionX(Form("dt_bg_rec%i",i),esub, esub+gsize); 

  int j=1;
  while(rbin<=1000){
  
    rbin = j*50;
    TH1D *gadd  = dtsingles->ProjectionX(Form("dt_add_rec%i_rebin%i",i,rbin),eadd, eadd+gsize); 
    gadd->Add(gsub,-1);
    gadd->Rebin(rbin);
    gadd->GetXaxis()->SetRangeUser(0,time);
    gadd->SetTitle(Form("150keV, first %.1fms, binwidth = %ius",time, rbin));
    c1->cd(j);
    gadd->Draw();
    
    TH1D *gadd1  = dtsingles->ProjectionX(Form("dt1_add_rec%i_rebin%i",i,rbin),eadd, eadd+gsize); 
    gadd1->Add(gsub,-1);
    gadd1->Rebin(rbin);
    gadd1->GetXaxis()->SetRangeUser(0.002*rbin,time);
    gadd1->SetTitle(Form("150keV, first %.1fms(no 1st bin), binwidth = %ius",time, rbin));
    c2->cd(j);
    gadd1->Draw();

    TH1D *gbg   = dtsingles->ProjectionX(Form("dt_bg_rec%i_rebin%i",i,rbin),eadd, eadd+gsize);
    gbg->Add(gsub,-1);
    gbg->Rebin(rbin);
    gbg->GetXaxis()->SetRangeUser(60,100);
    TF1 *fbg = new TF1("fbg","[0]");
    gbg->Fit(fbg); 
    
    TF1 *fun = new TF1("fun", "[0]*exp(-0.693/[1]*x)+[2]");
    fun->SetParameters(gadd->GetMaximum(),6,0.01);
    fun->FixParameter(2,fbg->GetParameter(0));
    gadd->Fit(fun);
    TText *chi2 = new TText();
    chi2 -> SetNDC();
    chi2 -> SetTextFont(1);
    chi2 -> SetTextColor(1);
    chi2 -> SetTextSize(0.05);
    chi2 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl = new TText();
    hl -> SetNDC();
    hl -> SetTextFont(1);
    hl -> SetTextColor(1);
    hl -> SetTextSize(0.05);
    hl -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));      
    gr->SetPoint(j-1,j,fun->GetParameter(1));
    gr->SetPointError(j-1,0,fun->GetParError(1));
    gr1->SetPoint(j-1,j,fun->GetChisquare()/fun->GetNDF());

    fun->SetParameters(gadd1->GetMaximum(),6,0.01);
    fun->FixParameter(2,fbg->GetParameter(0));
    gadd1->Fit(fun);
    TText *chi21 = new TText();
    chi21 -> SetNDC();
    chi21 -> SetTextFont(1);
    chi21 -> SetTextColor(1);
    chi21 -> SetTextSize(0.05);
    chi21 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl1 = new TText();
    hl1 -> SetNDC();
    hl1 -> SetTextFont(1);
    hl1 -> SetTextColor(1);
    hl1 -> SetTextSize(0.05);
    hl1 -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));
    gr2->SetPoint(j-1,j,fun->GetParameter(1));
    gr2->SetPointError(j-1,0,fun->GetParError(1));
    gr21->SetPoint(j-1,j,fun->GetChisquare()/fun->GetNDF());

    j++;

  }


  new TCanvas;
  gr->Draw("A*");
  new TCanvas;
  gr1->Draw("A*");
  new TCanvas;
  gr2->Draw("A*");
  new TCanvas;
  gr21->Draw("A*");

}


void NaDifBinWidth() {

  // bin count;
  double eadd = 882;
  double esub = 900;
  double gsize = 6;
  double time = 100; //unit:ms
  double starttime = 0.01; // unit:ms
  double startbin = 1*1000; // unit:us
  double tbin = time*1000;  // unit:us
  int rbin = 10;

  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdUnCor_slide535.root");
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  
  TCanvas *c1 = new TCanvas;
  c1->Divide(4,5);
  TCanvas *c2 = new TCanvas;
  c2->Divide(4,5);

  TGraphErrors *gr = new TGraphErrors(20);
  gr->SetTitle(Form("885keV, first %.1fms",time));
  TGraphErrors *gr1 = new TGraphErrors(20);
  gr1->SetTitle(Form("Chi2 885keV, first %.1fms",time));

  TGraphErrors *gr2 = new TGraphErrors(20);
  gr2->SetTitle(Form("885keV, first %.1fms(no 1st bin)",time));
  TGraphErrors *gr21 = new TGraphErrors(20);
  gr21->SetTitle(Form("Chi2 885keV, first %.1fms(no 1st bin)",time));
  
  int i = 5;
  TH2D *dtsingles = (TH2D *)f->Get(Form("dt_singles_recna%i",i));
  TH1D *dt   = dtsingles->ProjectionX(Form("dt_recna%i",i),eadd, eadd+gsize); 
  TH1D *gsub = dtsingles->ProjectionX(Form("dt_bg_recna%i",i),esub, esub+gsize); 

  int j=1;
  while(rbin<=1000){
  
    rbin = j*50;
    TH1D *gadd  = dtsingles->ProjectionX(Form("dt_add_rec%i_rebin%i",i,rbin),eadd, eadd+gsize); 
    gadd->Add(gsub,-1);
    gadd->Rebin(rbin);
    gadd->GetXaxis()->SetRangeUser(0,time);
    gadd->SetTitle(Form("150keV, first %.1fms, binwidth = %ius",time, rbin));
    c1->cd(j);
    gadd->Draw();
    
    TH1D *gadd1  = dtsingles->ProjectionX(Form("dt1_add_rec%i_rebin%i",i,rbin),eadd, eadd+gsize); 
    gadd1->Add(gsub,-1);
    gadd1->Rebin(rbin);
    gadd1->GetXaxis()->SetRangeUser(0.002*rbin,time);
    gadd1->SetTitle(Form("150keV, first %.1fms(no 1st bin), binwidth = %ius",time, rbin));
    c2->cd(j);
    gadd1->Draw();
    
    TF1 *fun = new TF1("fun", "[0]*exp(-0.693/[1]*x)+[2]");
    fun->SetParameters(gadd->GetMaximum(),8,0.01);
    fun->SetParLimits(2,0,10);
    gadd->Fit(fun);
    TText *chi2 = new TText();
    chi2 -> SetNDC();
    chi2 -> SetTextFont(1);
    chi2 -> SetTextColor(1);
    chi2 -> SetTextSize(0.05);
    chi2 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl = new TText();
    hl -> SetNDC();
    hl -> SetTextFont(1);
    hl -> SetTextColor(1);
    hl -> SetTextSize(0.05);
    hl -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));      
    gr->SetPoint(j-1,j,fun->GetParameter(1));
    gr->SetPointError(j-1,0,fun->GetParError(1));
    gr1->SetPoint(j-1,j,fun->GetChisquare()/fun->GetNDF());

    fun->SetParameters(gadd1->GetMaximum(),8,0.01);
    fun->SetParLimits(2,0,10);
    gadd1->Fit(fun);
    TText *chi21 = new TText();
    chi21 -> SetNDC();
    chi21 -> SetTextFont(1);
    chi21 -> SetTextColor(1);
    chi21 -> SetTextSize(0.05);
    chi21 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl1 = new TText();
    hl1 -> SetNDC();
    hl1 -> SetTextFont(1);
    hl1 -> SetTextColor(1);
    hl1 -> SetTextSize(0.05);
    hl1 -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));
    gr2->SetPoint(j-1,j,fun->GetParameter(1));
    gr2->SetPointError(j-1,0,fun->GetParError(1));
    gr21->SetPoint(j-1,j,fun->GetChisquare()/fun->GetNDF());

    j++;

  }


  new TCanvas;
  gr->Draw("A*");
  new TCanvas;
  gr1->Draw("A*");
  new TCanvas;
  gr2->Draw("A*");
  new TCanvas;
  gr21->Draw("A*");

}


void dt_gateTOF(){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/TOF_gate150_exdTOF.root");
  TH2D *dtTOF   = (TH2D *)f->Get("dt_TOF_30ms_150");
  TH2D *dtTOFbg = (TH2D *)f->Get("dt_TOF_30ms_BGR150");
  
  double binwindow = 70; //binwidth = 10;
  double bintof = 477;
  int rbin = 50;// orginal binwidth = 10us;
  
  TCanvas *c1 = new TCanvas;
  c1->Divide(4,4);
  TCanvas *c2 = new TCanvas;
  c2->Divide(4,4);
  TCanvas *c3 = new TCanvas;
  c3->Divide(4,4);
  
  TGraphErrors *gr = new TGraphErrors(16);
  gr->SetTitle("150keV, first 30ms");
  TGraphErrors *gr1 = new TGraphErrors(16);
  gr1->SetTitle("Chi2 150keV, first 30ms");

  TGraphErrors *grc = new TGraphErrors(16);
  grc->SetTitle("150keV, first 30ms, SubBG");
  TGraphErrors *grc1 = new TGraphErrors(16);
  grc1->SetTitle("Chi2 150keV, first 30ms, SubBG");
  for(int i=0;i<16;i++){
    bintof = 477 + i*23.;
    TH1D *dt   = dtTOF->ProjectionY(Form("dt_rec%i",i),bintof,bintof+binwindow);
    TH1D *cdt  = dtTOF->ProjectionY(Form("dt_clone_rec%i",i),bintof,bintof+binwindow);
    TH1D *dtbg = dtTOFbg->ProjectionY(Form("dt_bg_rec%i",i),bintof,bintof+binwindow);
    
    c1->cd(i+1);
    dt->Draw();
    c2->cd(i+1);
    dtbg->Draw();
    cdt->Add(dtbg,-1);
    c3->cd(i+1);
    cdt->Draw();    
    dt->Rebin(rbin);
    cdt->Rebin(rbin);
    dtbg->Rebin(rbin);

    TF1 *fun = new TF1("fun","[0]*exp(-0.693/[1]*x)+[2]");
    fun->SetParameters(cdt->GetMaximum(),5,0.1);
    fun->SetParLimits(2,0,1000);
    dt->Fit(fun);
    TText *chi2 = new TText();
    chi2 -> SetNDC();
    chi2 -> SetTextFont(1);
    chi2 -> SetTextColor(1);
    chi2 -> SetTextSize(0.05);
    chi2 -> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl = new TText();
    hl -> SetNDC();
    hl -> SetTextFont(1);
    hl -> SetTextColor(1);
    hl -> SetTextSize(0.05);
    hl -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));
    gr->SetPoint(i,i+1,fun->GetParameter(1));
    gr->SetPointError(i,0,fun->GetParError(1));
    gr1->SetPoint(i,i+1,fun->GetChisquare()/fun->GetNDF());

    cdt->Fit(fun);
    TText *chi21 = new TText();
    chi21-> SetNDC();
    chi21-> SetTextFont(1);
    chi21-> SetTextColor(1);
    chi21-> SetTextSize(0.05);
    chi21-> DrawText(0.5, 0.5, Form("X2 = %.2f",fun->GetChisquare()/fun->GetNDF()));
    TText *hl1 = new TText();
    hl1 -> SetNDC();
    hl1 -> SetTextFont(1);
    hl1 -> SetTextColor(1);
    hl1 -> SetTextSize(0.05);
    hl1 -> DrawText(0.5, 0.4, Form("hl = %.2f(%.3f)",fun->GetParameter(1), fun->GetParError(1)));
    grc->SetPoint(i,i+1,fun->GetParameter(1));
    grc->SetPointError(i,0,fun->GetParError(1));
    grc1->SetPoint(i,i+1,fun->GetChisquare()/fun->GetNDF());
  }

  new TCanvas;
  gr->Draw("A*");
  new TCanvas;
  gr1->Draw("A*");
  new TCanvas;
  grc->Draw("A*");
  new TCanvas;
  grc1->Draw("A*");


}


void TOF_gateDT() {
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/TOF_gate150_exdTOF.root");
  TH2D *dtTOF   = (TH2D *)f->Get("dt_TOF_30ms_150");
  TH2D *dtTOFbg = (TH2D *)f->Get("dt_TOF_30ms_BGR150");

  TCanvas *c1 = new TCanvas;
  c1->Divide(5,2);
  TCanvas *c2 = new TCanvas;
  c2->Divide(5,2);
  TCanvas *c3 = new TCanvas;
  c3->Divide(5,2);
  double startTbin = 1; // binwidth = 10us;
  double twindow = 3.*100; //time window = 3ms;
  int count = 1;
  while(startTbin<=30000) {
    TH1D *tof   = dtTOF->ProjectionX(Form("tof_at_%.1fms",startTbin/100),startTbin, startTbin+twindow-1);
    TH1D *ctof  = dtTOF->ProjectionX(Form("tof_clone_at_%.1fms",startTbin/100),startTbin, startTbin+twindow-1);
    TH1D *bgtof = dtTOFbg->ProjectionX(Form("tofBG_at_%.1fms",startTbin/100),startTbin, startTbin+twindow);
  
    c1->cd(count);
    tof->Draw();  
    c2->cd(count);
    bgtof->Draw();  
    ctof->Add(bgtof,-1);
    c3->cd(count);
    ctof->Draw();  
    tof->Draw("same");
    ctof->SetLineColor(kRed);
    tof  ->GetXaxis()->SetRangeUser(9000,15600);
    bgtof->GetXaxis()->SetRangeUser(9000,15600);
    ctof ->GetXaxis()->SetRangeUser(9000,15600);

    startTbin += twindow;
    count++; 
  }
  
  startTbin = 1;
  double twindow1 = 3.*100;
  TGraphErrors *grtof = new TGraphErrors(5);
  new TCanvas;
  for(int i=0;i<=100;i+=4){
    if(startTbin>2000) break;
    TH1D *tof   = dtTOF->ProjectionX(Form("Ctof_at_%.1fms",startTbin/100),startTbin, startTbin+twindow-1);
    TH1D *bgtof = dtTOFbg->ProjectionX(Form("CtofBG_at_%.1fms",startTbin/100),startTbin, startTbin+twindow);
    
    tof->Add(bgtof,-1);
    tof->Rebin(4);
    if(i==0) tof->Draw();  
    else{tof->SetLineColor(i/4); tof->Draw("same");}
    tof->GetXaxis()->SetRangeUser(9000,16000);
    
    //TH1D *ctof = (TH1D *)tof->Clone(Form("%s_clone",tof->GetName()));
    //ctof->Fit("gaus");
    //TF1 *g = (TF1*)ctof->GetListOfFunctions()->FindObject("gaus");
    //grtof->SetPoint(i,i/4+1,g->GetParameter(1));
    //grtof->SetPointError(i,0,g->GetParError(1));

    startTbin =1+ 100.*i;
  }

  TH2D *dt_TOF_885 = (TH2D *)f->Get("dt_TOF_60ms_885");
  TH2D *dt_TOF_BGR885 = (TH2D *)f->Get("dt_TOF_60ms_BGR885");
  int rbin = 10;
  TCanvas *c = new TCanvas;
  c->Divide(5,2);
  for(int i=0;i<10;i++){
    double tof = 12500.+i*230;
    tof=(tof-6000)/10.;
    TH1D *dt885 = dt_TOF_885->ProjectionY(Form("dt885_rec%i",i),tof,tof+70); 
    TH1D *dtbg885 = dt_TOF_BGR885->ProjectionY(Form("dtbg885_rec%i",i),tof,tof+70);
    TH1D *cdt885 = (TH1D *)dt885->Clone(Form("%s_clone",dt885->GetName()));
    cdt885->Add(dtbg885,-1);
    cdt885->Rebin(rbin);
    c->cd(i+1);
    cdt885->Draw();
    TF1 *fun = new TF1("fun","[0]*exp(-0.693/[1]*x)+[2]");
    fun->SetParameters(cdt885->GetMaximum(),9,0.01);
    cdt885->Fit(fun);
  }

}






