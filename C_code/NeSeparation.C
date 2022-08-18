



void NeSeparationMov(){

  std::ofstream ofile;
  ofile.open("NeSeparationMov1.txt");
  ofile << "Num/I:"
    <<setw(10)<< "PID_Int/D:"
    <<setw(10)<< "I2Pos_Int/D:"
    <<setw(10)<< "C150keV/D:"
    <<setw(10)<< "T150keV/D:"
    <<setw(10)<< "TError/D:" <<std::endl;


  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  TH2D *dt_single;
  TH2D *pid;
  TH2D *i2pos;
  TF1 *fx = new TF1("fx","[0]*exp(-0.693/[1]*x)+[2]",0,50);
  double t = 0; // half-life gated at 150keV;
  double terr = 0; // half-life error gated at 150keV;
  for(int row=i;row<16;i++){
    pid = (TH2D *)f->Get(Form("PID_2DMov1_%i",i));
    if(!pid){ //histogram doesn't exist
      ofile << i
        <<setw(10)<< 0.0
        <<setw(10)<< 0.0
        <<setw(10)<< 0.0
        <<setw(10)<< 0.0 
        <<setw(10)<< 0.0 <<std::endl;
      continue;
    } 
    double Ipid = pid->Integral(1,250,1,200);
    i2pos = (TH2D *)f->Get(Form("I2_TOF_2DMov1_%i",i));
    double Ii2pos = i2pos->Integral(1,190,1,300);

    dt_single = (TH2D *)f->Get(Form("dt_singleMov1_%i",i));
    TH1D *single = dt_single->ProjectionY(Form("single_%i",i),1,500);
    double count = single->Integral(147,152);//6 bins
    double bg = single->Integral(144,146); // 3 bin
    bg += single->Integral(153,155); // 3 bin
    count -= bg;
    if(count<5){
      t = 0.1;
      terr = 0.0;
    }else{
      TH1D *dt = dt_single->ProjectionX(Form("dt_%i",i),147,152);      
      dt->Rebin(50); // binwidth = 2ms;
      dt->GetXaxis()->SetRangeUser(0,50);
      fx->SetParameter(0,dt->GetBinContent(dt->GetMaximumBin()));   
      fx->SetParameter(1,5);
      fx->SetParameter(2,0.1);
      dt->Fit(fx);
      t = fx->GetParameter(1);
      terr = fx->GetParError(1);
    }
    ofile << i
      <<setw(10)<< Ipid
      <<setw(10)<< Ii2pos
      <<setw(10)<< count
      <<setw(10)<< t
      <<setw(10)<< terr <<std::endl;
  }




  return;
}

void NeSeparationRec(){

  std::ofstream ofile;
  ofile.open("NeSeparation0505.txt");
  ofile << "RowNum/I:"
    <<setw(10)<< "ColNum/I:"
    <<setw(10)<< "PID_Int/D:"
    <<setw(10)<< "I2Pos_Int/D:"
    <<setw(10)<< "C150keV/D:"
    <<setw(10)<< "T150keV/D:"
    <<setw(10)<< "TError/D:" <<std::endl;


  //TFile *f = TFile::Open("/home/zhu/packages/BCSSort/root_file/beta_NeSeparation/Ne2DCut1_Spectrum.root");
  TFile *f = TFile::Open("/home/zhu/packages/BCSSort/master_beta_prompt_op.root");
  TH2D *dt_single;
  TH2D *pid;
  TH2D *i2pos;
  TF1 *fx = new TF1("fx","[0]*exp(-0.693/[1]*x)+[2]",0,50);
  double t = 0; // half-life gated at 150keV;
  double terr = 0; // half-life error gated at 150keV;
  for(int row=0;row<5;row++){
    for(int col=0;col<5;col++){
      pid = (TH2D *)f->Get(Form("PID_2DCut3_%02i%02i",row,col));
      if(!pid){ //histogram doesn't exist
        ofile << row
          <<setw(10)<< col
          <<setw(10)<< 0.0
          <<setw(10)<< 0.0
          <<setw(10)<< 0.0
          <<setw(10)<< 0.0 
          <<setw(10)<< 0.0 <<std::endl;
        continue;
      } 
      double Ipid = pid->Integral(1,250,1,200);
      i2pos = (TH2D *)f->Get(Form("I2_TOF_2DCut3_%02i%02i",row,col));
      double Ii2pos = i2pos->Integral(1,190,1,300);

      dt_single = (TH2D *)f->Get(Form("dt_singleRec3_%02i%02i",row,col));
      TH1D *single = dt_single->ProjectionY(Form("single_%02i%02i",row,col),1,500);
      double count = single->Integral(147,152);//6 bins
      double bg = single->Integral(144,146); // 3 bin
      bg += single->Integral(153,155); // 3 bin
      count -= bg;
      if(count<5){
        t = 0.1;
        terr = 0.0;
      }else{
        TH1D *dt = dt_single->ProjectionX(Form("dt_%02i%02i",row,col),147,152);      
        dt->Rebin(50); // binwidth = 2ms;
        dt->GetXaxis()->SetRangeUser(0,50);
        fx->SetParameter(0,dt->GetBinContent(dt->GetMaximumBin()));   
        fx->SetParameter(1,5);
        fx->SetParameter(2,0.1);
        dt->Fit(fx);
        t = fx->GetParameter(1);
        terr = fx->GetParError(1);
      }
      ofile << row
        <<setw(10)<< col
        <<setw(10)<< Ipid
        <<setw(10)<< Ii2pos
        <<setw(10)<< count
        <<setw(10)<< t
        <<setw(10)<< terr <<std::endl;
    }
  }




  return;
}
