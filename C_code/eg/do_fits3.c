

void do_fits3(int rebin=10) {

  TH2* pw = (TH2*)_file0->Get("dt_TOF_60ms_885");
  TH2* bg = (TH2*)_file0->Get("dt_TOF_60ms_BGR885");

  pw = (TH2*)pw->Clone();
  bg = (TH2*)bg->Clone();

  //pw->RebinX(10);
  //bg->RebinX(10);
  pw->RebinY(rebin);
  bg->RebinY(rebin);

  int low,high;
  double width = 300;
  double shift = 100;
  double current = 12500;

  TH1 *p[16];
  TH1 *b[16];
  TH1 *pb[16];

  TF1 *f[16];

  gStyle->SetOptStat("nei");

  for(int x=0;x<16;x++) {
    low  = pw->GetXaxis()->FindBin(current+1);
    high = pw->GetXaxis()->FindBin(current+width-1);
    p[x] = pw->ProjectionY(Form("p%i",x),low,high); // from 11500 to 12000
    //p[x]->Rebin(10);
    b[x] = bg->ProjectionY(Form("b%i",x),low,high); // from 11500 to 12000
    //b[x]->Rebin(10);
    pb[x] = (TH1*)p[x]->Clone(Form("pb%i",x));
    pb[x]->Add(b[x],-1);
    pb[x]->SetTitle(Form("%i to %i",(int)current,int(current+width)));
    current+=shift;
  }





/*
  TF1 *f1  = new TF1("f1","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
  TF1 *f2  = new TF1("f2","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
  TF1 *fb1 = new TF1("fb1","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
  TF1 *fb2 = new TF1("fb2","[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);

  f1->SetParameters(25,5,5);
  f2->SetParameters(25,5,5);
  fb1->SetParameters(25,5,5);
  fb2->SetParameters(25,5,5);
 

  f1->SetParLimits(0,5,50);
  f1->SetParLimits(1,0,20);
  f1->SetParLimits(2,0,20);

  f2->SetParLimits(0,5,150);
  f2->SetParLimits(1,0,20);
  f2->SetParLimits(2,0,20);
  
  fb1->SetParLimits(0,5,50);
  fb1->SetParLimits(1,0,20);
  fb1->SetParLimits(2,0,20);

  fb2->SetParLimits(0,5,150);
  fb2->SetParLimits(1,0,20);
  fb2->SetParLimits(2,0,20);
*/

  pw->GetXaxis()->SetRangeUser(12500,15000);
  pw->GetXaxis()->SetLabelSize(0.08);
  pw->GetYaxis()->SetLabelSize(0.08);
  pw->SetStats(0);


  current = 12500;

  TCanvas *c = new TCanvas;
  c->SetTitle(Form("Rebin = %i",rebin));
  c->Divide(4,4);
  for(int x=0;x<16;x++) {
    
    f[x]  = new TF1(Form("f%i",x),"[0]*exp(-(0.69314718/[1])*x) + [2]",0,30);
    f[x]->SetParameters(25,5,5);
    f[x]->SetParLimits(0,5,2000);
    f[x]->SetParLimits(1,0,20);
    f[x]->SetParLimits(2,0,20);



    c->cd(x+1);
    pb[x]->Draw();
    pb[x]->SetStats(0);
    pb[x]->Fit(f[x]);

     c->Update();
     double y = gPad->GetUymax();

    TLatex l;
    l.SetTextSize(0.09);
    l.DrawLatex(15,y*.4,Form("t_{1/2} = %f",f[x]->GetParameter(1)));

    TPad *p = new TPad(Form("p%i",x),"",.6,.6,1,1);
    p->Draw();
    p->cd();
    pw->Draw("colz same");   

    c->Update();
 
    TLine *l1=new TLine(current,p->GetUymin(),current,p->GetUymax());
    TLine *l2=new TLine(current+width,p->GetUymin(),current+width,p->GetUymax());
    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);
    l1->Draw();
    l2->Draw();
    current+=shift;
 
  }


}
