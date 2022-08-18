


TList *cutlist(){
  TList *glist = new TList();
  TCutG *cutg;

double lna = 300;
  double wna = 1000;
  double horna = 100;
  double verna = 20;
  for(int i=0;i<18;i++){
    double start[2] = {12500,6200};
    start[0] += i*horna*1;
    start[1] -= i*verna*1;
    cutg = new TCutG(Form("recna%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]+wna);
    cutg->SetPoint(3,start[0],start[1]+wna);
    cutg->SetPoint(4,start[0],start[1]);
    glist->Add(cutg);
  }

  for(int i=0;i<22;i++){
    double start[2] = {10770,6140};
    start[0] += i*horna;
    start[1] -= i*verna;
    cutg = new TCutG(Form("rec%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+lna,start[1]);
    cutg->SetPoint(2,start[0]+lna,start[1]-wna);
    cutg->SetPoint(3,start[0],start[1]-wna);
    cutg->SetPoint(4,start[0],start[1]);
    glist->Add(cutg);
  }
  
  return glist;

}



void NadtDraw(int rbin=50){

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdUnCor_slide554NaOnly_dtbin1us.root");
  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  TH2 *PID = (TH2*)f->Get("pid_before"); 
  PID->SetTitle("PID Na chain"); 
  double eadd = 882;
  double esub = 900;
  double gsize = 5;
  double starttime = 0;
  double time = 100.;
  double tbin = time*100/rbin;
  double startbin = starttime*100/rbin;
  int tof = 12500;
  int width = 300;

  PID->SetStats(0);
  PID->GetXaxis()->SetRangeUser(12500,14520);
  PID->GetYaxis()->SetRangeUser(6000,7150);

  TCanvas *c1 = new TCanvas;
  c1->SetTitle(Form("rebin = %i",rbin));
  c1->Divide(4,3);
  TGraphErrors *gr = new TGraphErrors(12);
  gr->SetTitle(Form("half-life gated on 885keV, first %.1fms",time));
  gr->GetXaxis()->SetTitle("TOF at center of PID cut");
  gr->GetYaxis()->SetTitle("half-life(ms)");

  gStyle->SetOptStat(0);
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.5);
  
  TList *glist = cutlist();
  TCutG *cutg[18];

  
  for(int i=4;i<15;i++){
    tof = 12500 + i*100;    
    TH2 *pid = (TH2*)f->Get(Form("PID_recna%i",i));

    TH2 *dtsingles = (TH2D *)f->Get(Form("dt_singles_recna%i",i));
    dtsingles->RebinX(rbin);
    TH1 *singles = dtsingles->ProjectionY(Form("singles_recna%i",i),startbin,tbin);
    TH1 *dt = dtsingles->ProjectionX(Form("dt150_recna%i",i),eadd,eadd+gsize);
    TH1 *dtbg = dtsingles->ProjectionX(Form("dtbg_recna%i",i),esub, esub+gsize);

    if(i==12){
      c1->cd(1);
      singles->GetXaxis()->SetRangeUser(eadd-gsize, esub+2*gsize);
      TLine *l1 = new TLine(eadd,singles->GetMinimum(),eadd,singles->GetMaximum());
      TLine *l2 = new TLine(eadd+gsize,singles->GetMinimum(),eadd+gsize,singles->GetMaximum());
      TLine *l3 = new TLine(esub,singles->GetMinimum(),esub,singles->GetMaximum());
      TLine *l4 = new TLine(esub+gsize,singles->GetMinimum(),esub+gsize,singles->GetMaximum());
      l1->SetLineColor(kRed);
      l2->SetLineColor(kRed);
      l3->SetLineColor(kBlue);
      l4->SetLineColor(kBlue);
      singles->SetTitle("Gamma Energy Gate");
      singles->GetXaxis()->SetTitle("Energy(keV)");
      singles->GetYaxis()->SetTitle("Counts/keV");
      singles->Draw();
      l1->Draw("same");
      l2->Draw("same");
      l3->Draw("same");
      l4->Draw("same");
    }

    dt->SetTitle(Form("TOF = [%i, %i]",tof, tof+width));
    dt->GetXaxis()->SetTitle("Decay Time(ms)");
    dt->GetXaxis()->SetRangeUser(0,time);
    dt->GetYaxis()->SetTitle("Counts/ms");
    dt->GetXaxis()->SetTitleSize(0.05);
    dt->GetYaxis()->SetTitleSize(0.05);

    dt->Add(dtbg,-1);
    dt->GetXaxis()->SetRangeUser(0,time);
    c1->cd(i-2);
    dt->Draw();
    dt->SetStats(0);


    TF1 *fun = new TF1(Form("fun%i",i), "[0]*exp(-0.69314718/[1]*x)+[2]",0,100);
    //fun->SetParameters(dt->GetMaximum(),5,0.1);
    fun->SetParameters(50,8,0.1);
    fun->SetParLimits(0,5,2000);
    fun->SetParLimits(1,0,20);
    fun->SetParLimits(2,0,20);
    dt->Fit(fun);
    gr->SetPoint(i,tof+width/2,fun->GetParameter(1));
    gr->SetPointError(i,0,fun->GetParError(1));

    c1->Update();
    double y = gPad->GetUymax();

    TLatex l;
    l.SetTextSize(0.09);
    l.DrawLatex(25,y*.35,Form("t_{1/2} = %.2f(%.3f)ms",fun->GetParameter(1), fun->GetParError(1)));

    pid->SetStats(0);
    pid->GetXaxis()->SetRangeUser(12500,14520);
    pid->GetYaxis()->SetRangeUser(6000,7150);
    TPad *p = new TPad(Form("p%i",i),"",.6,.6,1,1);
    p->Draw();
    p->cd();
    PID->Draw("colz same");
    cutg[i] = (TCutG *)glist->FindObject(Form("recna%i",i));
    cutg[i]->SetLineColor(kRed);
    cutg[i]->SetLineWidth(2);
    cutg[i]->Draw("same");
    //pid->Draw("same"); 

  }
  new TCanvas;
  gr->GetXaxis()->SetRangeUser(12300,14150);
  gr->Draw("*A");


}



void NedtDraw(int flag=-1, double e=154){

  // bin count;
  double eadd = 147;
  //double esub = 155;
  double esub = e;
  double gsize = 5;
  double time = 30; //unit:ms
  double starttime = 0; // unit:ms
  double startbin = starttime*100; // unit:10us
  double tbin = time*100;  // unit:10us
  int rbin = 50;
  int tof = 10770;
  int width = 300;

  int size = 22;
  TF1 *fun[size];
  TH2 *pid[size];
  TH2 *dtsingles[size];
  TH1 *dt[size];
  TH1 *cdt[size];
  TH1 *dtbg[size];
  TH1 *singles[size];
  TList *glist = cutlist();
  TCutG *cutg[size];

  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/correlation1bestT/singles_TOFexdUnCor_slide515_dtbin1us.root");
  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  TH2 *PID = (TH2*)f->Get("pid_before");  
  PID->SetTitle("PID Ne chain"); 

  PID->SetStats(0);
  PID->GetXaxis()->SetRangeUser(10700,13400);
  PID->GetYaxis()->SetRangeUser(4700,6150);

  TCanvas *c1 = new TCanvas;
  c1->SetTitle("dt fit");
  c1->Divide(4,4);
  //TCanvas *c2 = new TCanvas;
  //c2->Divide(4,4);
  //TCanvas *c3 = new TCanvas;
  //c3->Divide(4,4);
  
  TGraphErrors *gr2 = new TGraphErrors(22);
  TGraphErrors *gr = new TGraphErrors(22);
  gr->SetTitle(Form("half-life gated on 150keV, first %.1fms",time));
  gr->GetXaxis()->SetTitle("TOF at center of PID cut");
  gr->GetYaxis()->SetTitle("half-life(ms)");

  gStyle->SetOptStat("nei");
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.5);


  for(int i=2;i<18;i++){
    tof = 10770 + i*100;
    pid[i] = (TH2*)f->Get(Form("PID_rec%i",i));

    dtsingles[i] = (TH2 *)f->Get(Form("dt_singles_rec%i",i));
    dtsingles[i]->RebinX(rbin);
    singles[i] = dtsingles[i]->ProjectionY(Form("singles_rec%i",i),startbin,tbin);
    dt[i] = dtsingles[i]->ProjectionX(Form("dt150_rec%i",i),eadd,eadd+gsize);
    cdt[i] = dtsingles[i]->ProjectionX(Form("cdt150_rec%i",i),eadd,eadd+gsize);
    dtbg[i] = dtsingles[i]->ProjectionX(Form("dtbg_rec%i",i),esub, esub+gsize*6);
    dtbg[i]->Scale(1/6.);
    //c2->cd(i-1);
    //cdt[i]->Draw();
    //c3->cd(i-1);
    //dtbg[i]->Draw("hist");

    dt[i]->SetTitle(Form("TOF = [%i, %i]",tof, tof+width));
    dt[i]->GetXaxis()->SetTitle("Decay Time(ms)");
    dt[i]->GetXaxis()->SetRangeUser(0,time);
    dt[i]->GetYaxis()->SetTitle("Counts/0.5ms");
    dt[i]->GetYaxis()->SetTitleSize(0.05);
    dt[i]->GetXaxis()->SetTitleSize(0.05);

    dt[i]->Add(dtbg[i],-1);
    dt[i]->GetXaxis()->SetRangeUser(0,time);
    c1->cd(i-1);
    dt[i]->Draw("hist");

    if(flag<0){
      fun[i] = new TF1(Form("fun%i",i), "[0]*exp(-0.69314718/[1]*x)+[2]",0,30);
      fun[i]->SetParameters(50,5,5);
      fun[i]->SetParLimits(0,0,200);
      fun[i]->SetParLimits(1,0,10);
      fun[i]->SetParLimits(2,0,50);
    }else{
      fun[i] = new TF1(Form("fun%i",i), "[0]*exp(-0.69314718/[1]*x)+[2]*exp(-0.69314718/[3]*x)+[4]");
      fun[i]->SetParameters(25,4,10,7,1);
      fun[i]->SetParLimits(0,0,2000);
      fun[i]->SetParLimits(1,2,5);
      fun[i]->SetParLimits(2,0,2000);
      fun[i]->SetParLimits(3,5,8.5);
      fun[i]->SetParLimits(4,0,200);
    }
    dt[i]->Fit(fun[i],"R");
    TF1 *fit = new TF1("fit","[0]*exp(-0.69314718/[1]*x)+[2]",0,30);
    fit->SetParameters(fun[i]->GetParameter(0),fun[i]->GetParameter(1),fun[i]->GetParameter(2));
    fit->SetLineColor(kRed);
    fit->Draw("same");
    gr->SetPoint(i,tof+width/2,fun[i]->GetParameter(1));
    gr->SetPointError(i,0,fun[i]->GetParError(1));
    if(flag>0){
      gr2->SetPoint(i,tof+width/2+100,fun[i]->GetParameter(3));
      gr2->SetPointError(i,0,fun[i]->GetParError(3));
    }

    c1->Update();
    double y = gPad->GetUymax();

    TLatex l;
    l.SetTextSize(0.1);
    //l.DrawLatex(15,y*.5,Form("A1=%.2f(%.3f) \t t_{1/2}=%.2f(%.3f)ms",fun[i]->GetParameter(0),fun[i]->GetParError(0),fun[i]->GetParameter(1), fun[i]->GetParError(1)));
    l.DrawLatex(13.5,y*.45,Form("t_{1/2} = %.2f(%.3f)ms",fun[i]->GetParameter(1), fun[i]->GetParError(1)));
    if(flag>0){
      TLatex ll;
      //ll.DrawLatex(13,y*.35,Form("t2_{1/2} = %.2f(%.3f)ms",fun[i]->GetParameter(3), fun[i]->GetParError(3)));
      ll.DrawLatex(13,y*.35,Form("A2=%.2f(%.3f) \t t_{1/2}=%.2f(%.3f)ms",fun[i]->GetParameter(2),fun[i]->GetParError(2),fun[i]->GetParameter(3), fun[i]->GetParError(3)));
    }

    pid[i]->SetStats(0);
    pid[i]->GetXaxis()->SetRangeUser(10700,13400);
    pid[i]->GetYaxis()->SetRangeUser(4700,6150);
    TPad *p = new TPad(Form("p%i",i),"",.6,.6,1,1);
    //p->SetLogz();
    p->Draw();
    p->cd();
    PID->Draw("colz same");
    cutg[i] = (TCutG *)glist->FindObject(Form("rec%i",i));
    cutg[i]->SetLineColor(kRed);
    cutg[i]->SetLineWidth(2);
    cutg[i]->Draw("same");

  }
  new TCanvas;
  gr->GetXaxis()->SetRangeUser(10700,13000);
  gr->Draw("A*");
  if(flag>0){
    new TCanvas;
    gr2->GetXaxis()->SetRangeUser(10700,13000);
    gr2->SetLineColor(kRed);
    gr2->Draw("A*");
  }

}




