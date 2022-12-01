

{

  //TH2D *dt_singles_ne30  = (TH2D *)_file0 ->Get("addback_ne30");
  //TH2D *dt_singles_ne31  = (TH2D *)_file0 ->Get("addback_ne31");
  TH2D *dt_singles_ne30  = (TH2D *)_file0 ->Get("tae7");
  TH2D *dt_singles_ne31  = (TH2D *)_file0 ->Get("tae6");

  TH2D *cdt0 = (TH2D *)dt_singles_ne30 ->Clone("ne30_adb" );
  TH2D *cdt1 = (TH2D *)dt_singles_ne31 ->Clone("ne31_adb" );

  int rebin = 8; // bin in 1ms;
  int time = 15;
  int tbins = time/100.*(1000/rebin);  // 20ms
  //int tbins = time/100.*(100/rebin);  // 20ms
  cdt0 ->RebinX(rebin);
  cdt1 ->RebinX(rebin);

  TH1D *ha0[150];
  TH1D *ha1[150];
  for(int i=0;i<tbins;i++){
    double tmin = i+1;
    double tmax = i+1;
    ha0[i]  = (TH1D *)cdt0 ->ProjectionY(Form("ne30_%i",i),tmin,tmax);
    ha1[i]  = (TH1D *)cdt1 ->ProjectionY(Form("ne31_%i",i),tmin,tmax);
  }

  double energy[16] = {150, 365, 373, 383, 423, 457, 602, 622, 698, 847, 1057, 1120, 1185, 2263, 3176, 3691};
  double lower[16]  = {147, 365, 373, 382, 420, 455, 600, 620, 693, 844, 1056, 1114, 1182, 2262, 3170, 3688};
  double upper[16]  = {152, 368, 375, 387, 427, 460, 605, 626, 701, 850, 1061, 1123, 1189, 2269, 3182, 3696};
  double bgmin[6]   = {153, 345, 377, 344, 378, 361};
  double bgmax[6]   = {158, 348, 380, 347, 381, 364};
 
  double     t[150];
  double    c0[150];
  double    c1[150];
  double   bg0[150]; 
  double   bg1[150]; 
  double  terr[150];
  double c0err[150];
  double c1err[150];
 
  TCanvas *ca1 = new TCanvas;
  TCanvas *ca2 = new TCanvas;
  ca1->Divide(5,4);  
  ca2->Divide(5,4);  

  int j=2;
  for(int i=0;i<tbins;i++){
    ca1->cd(i+1);ha0[i] ->GetXaxis()->SetRangeUser(lower[j]-20, upper[j]+20);ha0[i] ->Draw();    
    ca2->cd(i+1);ha1[i] ->GetXaxis()->SetRangeUser(lower[j]-20, upper[j]+20);ha1[i] ->Draw();    

    t[i]   = cdt0->GetXaxis()->GetBinCenter(i+1);
    c0[i]  = ha0[i] ->Integral(lower[j],upper[j]); 
    c1[i]  = ha1[i] ->Integral(lower[j],upper[j]); 
    
    //bg0[i] = ha0[i] ->Integral(bgmin[j],bgmax[j]);
    //bg1[i] = ha1[i] ->Integral(bgmin[j],bgmax[j]);
    
    c0[i] -= bg0[i]; 
    c1[i] -= bg1[i]; 
  
    terr[i] = 0;
    c0err[i]  = sqrt(c0[i])+sqrt(bg0[i]);
    c1err[i]  = sqrt(c1[i])+sqrt(bg1[i]);
    
    //cout << t[i] << "\t" << c0[i] << "\t" << terr[i] << "\t" << c0err[i] << std::endl;
    cout << t[i] << "\t" << c1[i] << "\t" << terr[i] << "\t" << c1err[i] << std::endl;

 
    ca1->Update();
    ca2->Update();
  } 
  //printf("hello\n");
  TGraphErrors *gr0  = new TGraphErrors(tbins, t,c0, terr, c0err);
  TGraphErrors *gr1  = new TGraphErrors(tbins, t,c1, terr, c1err);

  gr0  ->SetTitle(Form("ne30_%.1f", energy[j]));
  gr1  ->SetTitle(Form("ne31_%.1f", energy[j]));

  TCanvas *c = new TCanvas;
  c->Divide(1,2);
  c->cd(1);gr0 ->Draw("A*");
  c->cd(2);gr1 ->Draw("A*");

  TF1 *fx0  = new TF1("fx0", "[0]*exp(-0.693/[1]*x)+[2]");
  TF1 *fx1  = new TF1("fx1", "[0]*exp(-0.693/[1]*x)+[2]");
  fx0 ->SetParameters(20,2,-1);
  fx1 ->SetParameters(20,2,-1);
 
  fx0 ->SetParLimits(1,0,20);
  fx1 ->SetParLimits(1,0,20);
 
  gr0  -> Fit(fx0);
  gr1  -> Fit(fx1);

  std::cout<<"t1/2(30Ne) =  " <<fx0 ->GetParameter(1)<<std::endl;
  std::cout<<"t1/2(31Ne) =  " <<fx1 ->GetParameter(1)<<std::endl;

}
