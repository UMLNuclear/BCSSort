#include<map>


#include<TFile.h>
#include<TCutG.h>
#include<TH2D.h>
#include<TChain.h>
#include<TKey.h>
#include<TList.h>

#include<Implant.h>
#include<util.h>
#include<TChannel.h>
#include<DetHit.h>



void Process(std::vector<DetHit> vec){
  int sssd_flag = 0;
  std::vector<DetHit> FL;
  std::vector<DetHit> BL;
  DetHit pin1_i2s;
  DetHit pin1;
  DetHit i2s_i2n;
  double pin1_charge = -1;
  double pin1_cout = 0;
  for(size_t y=0;y<vec.size();y++){
    if(vec[y].GetNumber()==181){
      pin1_cout++;
      FillHistogram("PIN1E_in",10e3,0,10e3,vec[y].GetCharge());
      if(vec[y].GetCharge()>100) FillHistogram("PIN1E2_in",10e3,0,10e3,vec[y].GetCharge());
    }
    FillHistogram("PIN1_Count",10,0,10,pin1_cout);
    switch(vec[y].GetNumber()){
      case 40 ... 79:  //FL
        FL.push_back(vec[y]);
        break;

      case 120 ... 159:  //BL
        BL.push_back(vec[y]);
        break;

      case 160 ... 175:  //SSSD
        if(vec[y].GetCharge()>100){
          sssd_flag++;
        }
        break;

      case 177: //pin1_i2s
        pin1_i2s = vec[y];
        break;

      case 180: //i2s_i2n
        i2s_i2n = vec[y];
        break;

      case 181: //pin1
        pin1 = vec[y];
        pin1_charge = vec[y].GetCharge();
        break;

      default:
        break;
    }
  }
  if(pin1_charge>0){
    FillHistogram("PIN1E",10e3,0,10e3,pin1_charge);
    if(pin1_charge>100) FillHistogram("PIN1E2",10e3,0,10e3,pin1_charge);
    if(sssd_flag==0){
      if(FL.size()>0 && BL.size()>0){
        //FillHistogram("pid_isgood_nosssd",2e3,0,2e4,pin1_i2s.GetCharge(), 10e3,0,10e3,pin1.GetCharge());
        FillHistogram("PIN1E_isgood_nosssd", 10e3,0,10e3,pin1.GetCharge());
        if(pin1.GetCharge()>100) FillHistogram("PIN1E_isgood_nosssd2", 10e3,0,10e3,pin1.GetCharge());
      }
      FillHistogram("PIN1E_nosssd",10e3,0,10e3,pin1.GetCharge());
      if(pin1.GetCharge()>100) FillHistogram("PIN1E_nosss2",10e3,0,10e3,pin1.GetCharge());
    }  
  }


}




void listsort(){

  int faddress;
  int fnumber;
  double ftimestamp;
  double fcharge;
  gChain->SetBranchAddress("address",   &faddress);
  gChain->SetBranchAddress("number",    &fnumber);
  gChain->SetBranchAddress("timestamp", &ftimestamp);
  gChain->SetBranchAddress("charge",    &fcharge);
  TChannel::ReadDetMapFile();
  
  TList *gList1 = new TList;
  TList *gList2 = new TList;
  TList *gList3 = new TList;
  TList *gList4 = new TList;
  TList *gList5 = new TList;
  TList *gList6 = new TList;
  TH1 *hist1;
  TH1 *hist2;
  TH1 *hist3;
  TH1 *hist4;
  TH1 *hist5;
  TH1 *hist6;
  string hname1;
  string hname2;
  string hname3;
  string hname4;
  string hname5;
  string hname6;

  double starttime = 0;
  double buildtime = BUILDTIME;
  std::vector<DetHit> vec_hit;

  long n = gChain->GetEntries();
  //n = 1e6;
  long x = 0;
  for(x=0;x<n;x++){
    gChain->GetEntry(x);
    DetHit *fhit = new DetHit;
    fhit->SetAddress(faddress);
    fhit->SetNumber(fnumber);
    fhit->SetTimestamp(ftimestamp);
    fhit->SetCharge(fcharge);
    
    if(fhit->GetNumber()<208){ // LaBr
      int labr_num = fhit->GetNumber()- 192;
      hname1 = Form("Labr3_Eraw_%i", labr_num);
      hname2 = Form("Labr3_Ecal_%i", labr_num);
      hname3 = Form("Labr3_Time_%i", labr_num);
      hist1 = (TH1 *)gList1->FindObject(hname1.c_str());
      hist2 = (TH1 *)gList2->FindObject(hname2.c_str());
      hist3 = (TH1 *)gList3->FindObject(hname3.c_str());
      if(!hist1){
        hist1 = new TH1D(hname1.c_str(), hname1.c_str(), 65536,0,65536);
        gList1->Add(hist1);
      }
      if(!hist2){
        //hist2 = new TH1D(hname2.c_str(), hname2.c_str(), 4000,0,4000);
        hist2 = new TH1D(hname2.c_str(), hname2.c_str(), 3000,0,3000);
        gList2->Add(hist2);
      }
      if(!hist3){
        hist3 = new TH1D(hname3.c_str(), hname3.c_str(), 1440,0,86000);
        gList3->Add(hist3);
      }
      hist1->Fill(fhit->GetCharge());
      hist2->Fill(fhit->GetEnergy());
      hist3->Fill(fhit->GetTimestamp()/1e9);//unit:sec
    }else{ // HPGe
      int hpge_num = fhit->GetNumber() - 208;
      hname4 = Form("Clover_Eraw_%i", hpge_num);
      hname5 = Form("Clover_Ecal_%i", hpge_num);
      hname6 = Form("Clover_Time_%i", hpge_num);
      hist4 = (TH1 *)gList4->FindObject(hname4.c_str());
      hist5 = (TH1 *)gList5->FindObject(hname5.c_str());
      hist6 = (TH1 *)gList6->FindObject(hname6.c_str());
      if(!hist4){
        hist4 = new TH1D(hname4.c_str(), hname4.c_str(), 65536,0,65536);
        gList4->Add(hist4);
      }
      if(!hist5){
        hist5 = new TH1D(hname5.c_str(), hname5.c_str(), 3000,0,3000);
        //hist5 = new TH1D(hname5.c_str(), hname5.c_str(), 4000,0,4000);
        gList5->Add(hist5);
      }
      if(!hist6){
        hist6 = new TH1D(hname6.c_str(), hname6.c_str(), 1440,0,86000);
        gList6->Add(hist6);
      }
      hist5->Fill(fhit->GetEnergy());
      hist4->Fill(fhit->GetCharge());
      hist6->Fill(fhit->GetTimestamp()/1e9);//unit:sec
    }

    if((x%50000)==0){
      printf("  on entry %lu / %lu     \r",x,n);
      fflush(stdout);
    }
  }
  printf("    on entry %lu / %lu   \n",x,n);
  gList1->Sort();
  gList2->Sort();
  gList3->Sort();
  gList4->Sort();
  gList5->Sort();
  gList6->Sort();

  std::string runnum = GetRunNumber(gChain->GetCurrentFile()->GetName());
  TFile *f = new TFile(Form("run%s_spectrum.root",runnum.c_str()),"recreate");
  TDirectory *current = gDirectory;
  TDirectoryFile *Clovers = new TDirectoryFile("Clovers", "Clovers");
  TDirectoryFile *LaBr3s = new TDirectoryFile("Labr3s", "Labr3s");
  f->cd("Clovers");
  TDirectory *Raw = new TDirectoryFile("Raw", "Raw");
  TDirectory *Cal = new TDirectoryFile("Cal", "Cal");
  TDirectory *Time = new TDirectoryFile("Time", "Time");
  f->cd("Clovers/Raw");
  gList4->Write(); 
  f->cd("Clovers/Cal");
  gList5->Write(); 
  f->cd("Clovers/Time");
  gList6->Write(); 
  f->cd("Labr3s");
  TDirectory *Raw2  = new TDirectoryFile("Raw2", "Raw2");
  TDirectory *Cal2  = new TDirectoryFile("Cal2", "Cal2");
  TDirectory *Time2 = new TDirectoryFile("Time2","Time2");
  f->cd("Labr3s/Raw2");
  gList1->Write(); 
  f->cd("Labr3s/Cal2");
  gList2->Write(); 
  f->cd("Labr3s/Time2");
  gList3->Write(); 
  f->Close();

  return;
}
