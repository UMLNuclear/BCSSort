

double background(double *x, double *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double lorentzianPeak(double *x, double *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}

double totalFunction(double *x, double *par) {
  return background(x,par) + lorentzianPeak(x,&par[3]);
}


TF1 *f=0;
TF1 *l=0;
TF1 *bg=0;

void drawFunction(double p0, double p1,double p2,double p3, double p4, double p5) {
  if(!f) f = new TF1("f",totalFunction,0,100,6);
  f->SetParameters(p0,p1,p2,p3,p4,p5);

  if(!bg) bg = new TF1("bg",background,0,100,3);
  bg->SetParameters(p0,p1,p2);
  bg->SetLineColor(kBlue);

  if(!l) l = new TF1("l",lorentzianPeak,0,100,3);
  l->SetParameters(p3,p4,p5);
  l->SetLineColor(kGreen);


  f->Draw();  
  l->Draw("same");  
  bg->Draw("same");  


}

void some_junk() {
  drawFunction(20,5,-.1,10000,10,50);
}



