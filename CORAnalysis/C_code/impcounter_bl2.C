
#define SA 0


typedef std::pair<int, int> pixel;

TCutG *ban;
TCutG *cutstrip;
std::map<pixel, TCutG *>hgcut;

//============= veto for imp =============//
bool veto(DetHit lgf, std::vector<DetHit> sssd){

  int    lgfn = lgf.GetStrip()-1;
  double lgfe = lgf.GetEnergy();
  for(int i=0;i<sssd.size();i++){
    int sssdn  = sssd[i].GetStrip()-1;
    long sssdc = sssd[i].GetCharge();
    if(cutstrip->IsInside(sssdn, lgfn)){
      if(ban->IsInside(sssdc, lgfe)){
        return true;
      }
    }
  }
  return false;
}

//=========== Initialize Map ===========//
void InitMap(std::map<pixel, Implant*> *mp1, std::map<pixel, Implant*> *mp2, std::map<pixel, std::vector<Decay> *> *mp3){
  std::map<pixel, Implant*> temp1;
  std::map<pixel, Implant*> temp2;
  std::map<pixel, std::vector<Decay> *> temp3;
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      temp1[pix] = new Implant;
      temp2[pix] = new Implant;
      temp3[pix] = new std::vector<Decay>;
    }
  }
  *mp1 = temp1;
  *mp2 = temp2;
  *mp3 = temp3;
}

//=========== Clean Map ===============//
void CleanMap(std::map<pixel, Implant*> *mp, double pin1t, int flush=0, std::vector<pixel> *vecpix=0){
  for(auto &it: (*mp)){
    if(it.second->fPIN1T<0) continue;
    double dt = pin1t - it.second->fPIN1T;
    if(dt>1e9){
      if(flush){
        if(vecpix) vecpix->push_back(it.first);
      }else{
        it.second->Clear();
      }
    }
  }
}

//============ After loop, clean all implants on map ================//
void CleanMapAll(std::map<pixel, Implant*> *mp, std::vector<pixel> *vecpix){
  for(auto &it:(*mp)){
    if(it.second->fPIN1T<0) continue;
    vecpix->push_back(it.first);
  }
}

//============== check map and return info ===============//
int check(std::map<pixel, Implant *> mp, pixel pix, int x=0, std::vector<pixel> *vecpix=0){
  int c = 0;
  for(int m=-x;m<=x;m++){
    for(int n=-x;n<=x;n++){  
      pixel temp = std::make_pair(pix.first+m, pix.second+n);
      if(temp.first<0 || temp.first>39) continue;
      if(temp.second<0 || temp.second>39) continue;
      if(mp[temp]->fPIN1T>0) {
        c++;
        if(vecpix) vecpix->push_back(temp);
      }
    }
  }
  return c;
}

//=============== insert BCSEvent to map =================//
void Insert(std::map<pixel, Implant*> *mp, pixel pix, Implant *imphits){
  (*mp)[pix]->Clear();
  imphits->Copy(*((*mp)[pix]));
}

//============= DecayPixel ================//
pixel DecayPixel(std::vector<DetHit> hgf, std::vector<DetHit>hgb){
  std::multimap<double, DetHit> mmhgf;
  std::multimap<double, DetHit> mmhgb;
  for(int i=0;i<hgf.size();i++){
    double hgfc = hgf[i].GetCharge();
    if(hgfc<10 || hgfc>25000) continue;
    mmhgf.insert(std::make_pair(hgfc, hgf[i]));
  }
  for(int i=0;i<hgb.size();i++){
    double hgbc = hgb[i].GetCharge();
    if(hgbc<10 || hgbc>25000) continue;
    mmhgb.insert(std::make_pair(hgbc, hgb[i]));
  }
  std::multimap<double, DetHit>::reverse_iterator rit1;
  for(rit1=mmhgf.rbegin();rit1!=mmhgf.rend();++rit1){
    double hgft = rit1->second.GetTimestamp();
    double hgfc = rit1->second.GetCharge();
    int    hgfn = rit1->second.GetStrip()-1;
    std::multimap<double, DetHit>::reverse_iterator rit2;
    for(rit2=mmhgb.rbegin();rit2!=mmhgb.rend();++rit2){
      double hgbt = rit2->second.GetTimestamp();
      double hgbc = rit2->second.GetCharge();
      int    hgbn = rit2->second.GetStrip()-1;
      double dt = hgft-hgbt;
      if(dt>=30 || dt<140){
        pixel cutpix = std::make_pair(hgfn/16+1, hgbn/16+1);
        if(hgcut[cutpix]->IsInside(hgfc,hgbc)){
          return std::make_pair(hgfn,hgbn);
        }
      }
    }
  }
  return std::make_pair(-1,-1);
}

//========== Correlation =============//
pixel Correlation(std::map<pixel, Implant* >mimp, pixel pix, double dect, int x=0){
  for(int m=-x;m<=x;m++){
    for(int n=-x;n<=x;n++){
      pixel temp = std::make_pair(pix.first+m, pix.second+n);
      if(temp.first<0 || temp.first>39) continue;
      if(temp.second<0 || temp.second>39) continue;
      if(mimp[temp]->fPIN1T>0){
        double dt = dect - mimp[temp]->fPIN1T;
        if(dt<0 || dt>1e9) continue;
        return temp; 
      }
    }
  }
  return std::make_pair(-1,-1);
}


//================== "main" code for running ===============//
void counter(int runnum=1040) {
  TChannel::ReadDetMapFile();
  
  int count0 = 0;
  int count1 = 0;
  int count2 = 0;
  int count3 = 0;
 
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate2.root");
  ban = (TCutG *)cutf2->Get("ban");

  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  cutstrip = (TCutG *)cutf3->Get("cutstrip");
  
  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/datafile/event/implant/event%i*root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long n = chan->GetEntries();
  long x = 0;

  Implant *fimp = new Implant;

  std::map<pixel, Implant*>fImpMap;
  std::map<pixel, Implant*>fBLMap;
  std::map<pixel, std::vector<Decay> *> fDecMap;

  //Initialize map
  InitMap(&fImpMap, &fBLMap, &fDecMap);

  for(x=0;x<n;x++){
    fevent->Clear();
    chan->GetEntry(x);
    bool flag_veto = false;
    if(fevent->LGFSize()==0 || fevent->LGBSize()==0 || fevent->LGFMax().GetEnergy()<0 || fevent->LGBMax().GetEnergy()<0) continue;
    //veto check
    if(fevent->SSSDSize()>0){ flag_veto = veto(fevent->LGFMax(), fevent->SSSD()); }

    // out of the veto
    if(flag_veto) continue;
    fimp->Clear();
    fimp->Set(&(fevent->fHits));
    double pin1t = fimp->fPIN1T;
    pixel piximp = fimp->GetPixel();
    count0++;
    //Correlation
    
    //Clean Map
    std::vector<pixel> vecpiximp;
    CleanMap(&fImpMap, pin1t, 1, &vecpiximp);
    for(int i=0;i<vecpiximp.size();i++){
      pixel temp = vecpiximp[i];
      //printf("impentry = %lu \t decaysize = %i\n", ximp,fbeta->DecaySize());
      fImpMap[temp]->Clear();
      count1++;
    }
    CleanMap(&fBLMap, pin1t);

    if(check(fBLMap, piximp, SA)){//reset the bl clock
      Insert(&fBLMap, piximp, fimp);
      count2++;
    }else{
      std::vector<pixel> vecpix;
      int found = check(fImpMap, piximp, SA, &vecpix);
      if(found){ //back to back
        fImpMap[piximp]->Clear();
        Insert(&fBLMap, piximp, fimp);
        count3++;
        for(int i=0;i<found;i++){
          pixel temp = vecpix[i];
          if(temp.first==piximp.first && temp.second==piximp.second) {count3++;continue;}
          Insert(&fBLMap, temp, fImpMap[temp]);
          fImpMap[temp]->Clear();
          count3++;
        }
      }else{
        Insert(&fImpMap, piximp, fimp);
      }
    }
    // check time written correctly      
    for(auto it:fBLMap){
      if(it.second->fPIN1T<0) continue;
      if(it.second->fPIN1T < pin1t) continue;
      if(it.first.first == piximp.first && it.first.second == piximp.second) continue;
      printf("BLMap time written wrong!!!\n");
    } 
    for(auto it:fImpMap){
      if(it.second->fPIN1T<0) continue;
      if(it.second->fPIN1T < pin1t) continue;
      if(it.first.first == piximp.first && it.first.second == piximp.second) continue;
      printf("ImpMap time written wrong!!!\n");
    } 


    if((x%5000)==0){
      printf("run%i.root: on entry %lu / %lu \r", runnum,x,n);
      fflush(stdout);
    }
  }
 
  //to collect all decay events

 
  //clean all fImpMap
  std::vector<pixel> vecpiximp;
  CleanMapAll(&fImpMap, &vecpiximp);
  for(int i=0;i<vecpiximp.size();i++){
    pixel temp = vecpiximp[i];
    fImpMap[temp]->Clear();
    count1++;
  }
  printf("run%i.root: on entry %lu / %lu \n", runnum,x,n);
 
  std::cout << "good imp     = " << count0 <<std::endl;
  std::cout << "fill in beta = " << count1 <<std::endl;
  std::cout << "3rd block imp= " << count2 <<std::endl;
  std::cout << "back to back = " << count3 <<std::endl;

}

