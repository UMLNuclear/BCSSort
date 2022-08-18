
void *junk(){
//TList *junk(){

  TList *gList = new TList();
  TCutG *cutg;
  
  double l = 700;
  double w = 800;
  double hor = 230;
  double ver= 40;

  for(int i=0;i<16;i++){
    double start[2] = {11100,5300};
    start[0] += i*hor;
    start[1] -= i*ver;
    cutg = new TCutG(Form("rec%i",i),5);
    cutg->SetPoint(0,start[0],start[1]);
    cutg->SetPoint(1,start[0]+l,start[1]);
    cutg->SetPoint(2,start[0]+l,start[1]+w);
    cutg->SetPoint(3,start[0],start[1]+w);
    cutg->SetPoint(4,start[0],start[1]);
    //gList->Add(cutg);  
  }
  cutg = new TCutG("rec16",5); 
  cutg->SetPoint(0,11100,5300);
  cutg->SetPoint(1,13200,4940);
  cutg->SetPoint(2,13200,5740);
  cutg->SetPoint(3,11100,6100);
  cutg->SetPoint(4,11100,5300);
  //gList->Add(cutg);  

  //==== PID cuts for Na chain ====//
  double lna = 300;
  double wna = 1000;
  //double horna = 230;
  //double verna = 40;
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
    cutg->SetLineColor(i+1);
    cutg->SetLineWidth(4);
    //cutg->Draw("same");
    //gList->Add(cutg);
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
    cutg->SetLineColor(i+1);
    cutg->SetLineWidth(2);
    //if(i>1 && i<18)
    //cutg->Draw("same");
    //gList->Add(cutg);
  }
  //for(int i=9;i<17;i++){
  //  double start[2] = {10770,5340};
  //  start[0] -= (i-8)*horna;
  //  //start[1] -= i*verna;
  //  cutg = new TCutG(Form("rec%i",i),5);
  //  cutg->SetPoint(0,start[0],start[1]);
  //  cutg->SetPoint(1,start[0]+lna,start[1]);
  //  cutg->SetPoint(2,start[0]+lna,start[1]+wna);
  //  cutg->SetPoint(3,start[0],start[1]+wna);
  //  cutg->SetPoint(4,start[0],start[1]);
  //  cutg->SetLineColor(i);
  //  cutg->SetLineWidth(2);
  //  cutg->Draw("same");
  //  gList->Add(cutg);
  //}

  TCutG *NaCut = new TCutG("NaCut",5);
  NaCut->SetPoint(0,12200,6300);
  NaCut->SetPoint(1,12200,7300);
  NaCut->SetPoint(2,15600,6600);
  NaCut->SetPoint(3,15600,5600);
  NaCut->SetPoint(4,12200,6300);
  NaCut->SetLineColor(kBlack);
  NaCut->SetLineWidth(3);
  //NaCut->Draw("same");

  TCutG *NeCut = new TCutG("NeCut",5);
  NeCut->SetPoint(0,10970,5100);
  NeCut->SetPoint(1,11970,5000);
  NeCut->SetPoint(2,11970,5000+wna);
  NeCut->SetPoint(3,10970,5100+wna);
  NeCut->SetPoint(4,10970,5100);
  NeCut->SetLineColor(kRed);
  NeCut->SetLineWidth(3);
  NeCut->Draw("same");
  
  TCutG *NeCut0 = new TCutG("NeCut0",5);
  NeCut0->SetPoint(0,11970,5000);
  NeCut0->SetPoint(1,12870,4820);
  NeCut0->SetPoint(2,12870,4820+wna);
  NeCut0->SetPoint(3,11970,5000+wna);
  NeCut0->SetPoint(4,11970,5000);
  NeCut0->SetLineColor(kYellow);
  NeCut0->SetLineWidth(3);
  NeCut0->Draw("same");

  return gList;

}
