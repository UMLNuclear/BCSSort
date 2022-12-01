

{
  double energy[16] = {150, 305, 373, 383, 423, 457, 602, 622, 698, 847, 1057, 1120, 1185, 2263, 3176, 3691};
  double lower[16]  = {147, 304, 373, 382, 420, 455, 600, 620, 693, 844, 1056, 1114, 1182, 2262, 3170, 3688};
  double upper[16]  = {152, 308, 376, 387, 427, 460, 605, 626, 701, 850, 1061, 1123, 1189, 2269, 3182, 3696};
  
  TH2D *h2[7];
  TH2D *ch2[7];

  int rebin = 8; // bin in 1ms;
  int time = 25;  // fitting range = 15ms;
  int tbins = time/100.*(1000/rebin);
  
  TH1D  *h1[7][42];
  double     t[42];
  double    c0[42];
  double    c1[42];
  double    c2[42];
  double    c3[42];
  double    c4[42];
  double    c5[42];
  double    c6[42];
  double    te[42];
  double   c0e[42];
  double   c1e[42];
  double   c2e[42];
  double   c3e[42];
  double   c4e[42];
  double   c5e[42];
  double   c6e[42];
  
  for(int i=0;i<7;i++){
    //h2[i] = (TH2D *)_file0->Get(Form("dt30_addback_31ne_%i",i));
    h2[i] = (TH2D *)_file0->Get(Form("recae%i",i));
    ch2[i] = (TH2D *)h2[i]->Clone(Form("%s_c",h2[i]->GetName()));
    ch2[i] -> RebinX(rebin);
    for(int j=0;j<tbins;j++){
      h1[i][j] = (TH1D *)ch2[i]->ProjectionY(Form("31ne_gate%i_%i",i,j),j+1,j+1);
    }
  }


  for(int i=0;i<tbins;i++){
    t[i] = ch2[0]->GetXaxis()->GetBinCenter(i+1);
    te[i]= 0;
    c0[i] = h1[0][i]->Integral(lower[2],upper[2]);
    c1[i] = h1[1][i]->Integral(lower[2],upper[2]);
    c2[i] = h1[2][i]->Integral(lower[2],upper[2]);
    c3[i] = h1[3][i]->Integral(lower[2],upper[2]);
    c4[i] = h1[4][i]->Integral(lower[2],upper[2]);
    c5[i] = h1[5][i]->Integral(lower[2],upper[2]);
    c6[i] = h1[6][i]->Integral(lower[2],upper[2]);
    c0e[i] = sqrt(c0[i]);
    c1e[i] = sqrt(c1[i]);
    c2e[i] = sqrt(c2[i]);
    c3e[i] = sqrt(c3[i]);
    c4e[i] = sqrt(c4[i]);
    c5e[i] = sqrt(c5[i]);
    c6e[i] = sqrt(c6[i]);
    std::cout<< t[i] << "\t" << c0[i] << "\t" << te[i] << "\t" <<c2e[i] <<std::endl;
  } 
  
  TGraphErrors *gr[7];
  for(int i=0;i<7;i++){
    switch (i){
      case 0: gr[i] = new TGraphErrors(tbins,t,c0,te,c0e);break;
      case 1: gr[i] = new TGraphErrors(tbins,t,c1,te,c1e);break;
      case 2: gr[i] = new TGraphErrors(tbins,t,c2,te,c2e);break;
      case 3: gr[i] = new TGraphErrors(tbins,t,c3,te,c3e);break;
      case 4: gr[i] = new TGraphErrors(tbins,t,c4,te,c4e);break;
      case 5: gr[i] = new TGraphErrors(tbins,t,c5,te,c5e);break;
      case 6: gr[i] = new TGraphErrors(tbins,t,c6,te,c6e);break;
      default: break;
    }
  }

  TF1 *fx[7];
  for(int i=0;i<7;i++){
    gr[i] -> SetTitle(Form("Ne31_gate%i",i));
    fx[i] = new TF1(Form("fx%i",i),"[0]*exp(-0.693/[1]*x)+[2]");
    fx[i] -> SetParameters(15,3,2);
    fx[i] -> SetParLimits(1,0,10);
  }
  
  TCanvas *c = new TCanvas;
  c->Divide(4,2);
  for(int i=0;i<7;i++){
    c->cd(i+1);
    gr[i]->Draw("A*");
    gr[i]->Fit(fx[i]);

  }

  for(int i=0;i<7;i++){
    std::cout<<"t1/2(31ne_"<<i<<") = "<< fx[i]->GetParameter(1)<<std::endl;
  }

}
