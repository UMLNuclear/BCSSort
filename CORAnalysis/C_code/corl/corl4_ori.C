
#define SA 1


typedef std::pair<int, int> pixel;

TCutG *ban;
TCutG *cutstrip;
TCutG *blob[6];
std::map<pixel, TCutG *>hgcut;

int countt = 0;
int count0 = 0; // counts out of veto;
int count1 = 0; // good implant will be written in TTree;
int count2 = 0; // black list: back to back;
int count3 = 0; // black list: pixel is still cooling down;
int count4 = 0; // multiple imp correlate to one decay;

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
void InitMap(std::map<pixel, BCSEvent*> *mp1, std::map<pixel, BCSEvent*> *mp2, std::map<pixel, std::vector<Decay> *> *mp3){
  std::map<pixel, BCSEvent*> temp1;
  std::map<pixel, BCSEvent*> temp2;
  std::map<pixel, std::vector<Decay> *> temp3;
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      temp1[pix] = new BCSEvent;
      temp2[pix] = new BCSEvent;
      temp3[pix] = new std::vector<Decay>;
    }
  }
  *mp1 = temp1;
  *mp2 = temp2;
  *mp3 = temp3;
}

//=========== Clean Map ===============//
void CleanMap(std::map<pixel, BCSEvent*> *mp, double pin1t, int flush=0, std::map<pixel, std::vector<Decay> *> *mpdec=0){
  for(auto &it: (*mp)){
    if(it.second->Pin1T()<0) continue;
    double dt = pin1t - it.second->Pin1T();
    if(dt>1e9){
      if(flush){
        count1++;
      }
      it.second->Clear();
    }
  }
}

//============ After loop, clean all implants on map ================//
void CleanMapAll(std::map<pixel, BCSEvent*> *mp){
  for(auto &it:(*mp)){
    if(it.second->Pin1T()<0) continue;
    count1++;
    it.second->Clear();
  }
}

//============== check map and return info ===============//
int check(std::map<pixel, BCSEvent *> mp, pixel pix, int x=0, std::vector<pixel> *vecpix=0){
  int c = 0;
  for(int m=-x;m<=x;m++){
    for(int n=-x;n<=x;n++){  
      pixel temp = std::make_pair(pix.first+m, pix.second+n);
      if(temp.first<0 || temp.first>39) continue;
      if(temp.second<0 || temp.second>39) continue;
      if(mp[temp]->Pin1T()>0) {
        c++;
        if(vecpix) vecpix->push_back(temp);
      }
    }
  }
  return c;
}

//=============== insert BCSEvent to map =================//
void Insert(std::map<pixel, BCSEvent*> *mp, pixel pix, std::vector<DetHit> hits){
  (*mp)[pix]->Clear();
  (*mp)[pix]->fHits = hits;
}

//============= DecayPixel ================//
pixel DecayPixel(std::vector<DetHit> hgf, std::vector<DetHit>hgb){
  int indexf = -1;
  int indexb = -1;
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
  std::multimap<double, DetHit>::reverse_iterator rit2;
  for(rit1=mmhgf.rbegin();rit1!=mmhgf.rend();++rit1){
    double hgft = rit1->second.GetTimestamp();
    double hgfc = rit1->second.GetCharge();
    int    hgfn = rit1->second.GetStrip()-1;
    for(rit2=mmhgb.rbegin();rit2!=mmhgb.rend();++rit2){
      double hgbt = rit2->second.GetTimestamp();
      double hgbc = rit2->second.GetCharge();
      int    hgbn = rit2->second.GetStrip()-1;
      double dt = hgft-hgbt;
      if(dt>=30 || dt<140){
        pixel cutpix = std::make_pair(hgfn/16+1, hgbn/16+1);
        if(hgcut[cutpix]->IsInside(hgfc,hgbc)){
          return std::make_pair(indexf,indexb);
        }
      }
    }
  }
  return std::make_pair(indexf,indexb);
}

//========== Correlation =============//
pixel Correlation(std::map<pixel, BCSEvent * >mimp, pixel pix, double dect, int x=0){
  for(int m=-x;m<=x;m++){
    for(int n=-x;n<=x;n++){
      pixel temp = std::make_pair(pix.first+m, pix.second+n);
      if(temp.first<0 || temp.first>39) continue;
      if(temp.second<0 || temp.second>39) continue;
      if(mimp[temp]->Pin1T()>0){
        double dt = dect - mimp[temp]->Pin1T();
        if(dt<0 || dt>1e9) continue;
        c++; 
      }
    }
  }
  if(c>1) count4++;
  return std::make_pair(-1,-1);
}


//================== "main" code for running ===============//
void counter(int runnum = 1040){

  TChannel::ReadDetMapFile();
 
  std::vector<double> tofpar = TOFCorrection::Get()->ReadFile(runnum);
  std::vector<double> tofpar1= TOFCorrection::Get()->ReadFile(runnum,"/home/zhu/packages/BCSSort/config/TOF/TOF_beta_offset.txt");
  if(tofpar.empty()){
    printf("run%i.root doesn't have TOF correction parameters\n",runnum);
    return;
  }

  TFile *cutf1 = TFile::Open("/home/zhu/packages/BCSSort/root_file/cut/pid_cut_allblobs.root");
  for(int i=0;i<6;i++){
    blob[i] = (TCutG *)cutf1->Get(Form("blob%i",i));
  }

  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate2.root");
  ban = (TCutG *)cutf2->Get("ban");

  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  cutstrip = (TCutG *)cutf3->Get("cutstrip");

  TFile *cutf4 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/decay_hgcc_cut.root");
  for(int i=1;i<4;i++){
    for(int j=1;j<4;j++){
      pixel cutpix = std::make_pair(i,j);
      hgcut[cutpix] = (TCutG*)cutf4->Get(Form("cut%i%i",i,j));
    }
  }
  
  //int runnum = 1040;
  TChain *chimp = new TChain("event");
  TChain *chdec = new TChain("event");
  chimp->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/implant/event%i*root",runnum));
  chdec->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/decay/highgain/event%i*.root",runnum));
  BCSEvent *fimp = new BCSEvent;
  BCSEvent *fdec = new BCSEvent;
  chimp->SetBranchAddress("BCSEvent", &fimp);
  chdec->SetBranchAddress("BCSEvent", &fdec);
  long nimp = chimp->GetEntries();
  long ndec = chdec->GetEntries();
  long ximp = 0;
  long xdec = 0;

  Implant *ffimp = new Implant;
  Decay *ffdec = new Decay;

  std::map<pixel, BCSEvent*>fImpMap;
  std::map<pixel, BCSEvent*>fBLMap;
  std::map<pixel, std::vector<Decay> *> fDecMap;

  //Initialize map
  InitMap(&fImpMap, &fBLMap, &fDecMap);

  for(ximp=0;ximp<nimp;ximp++){
    fimp->Clear();
    chimp->GetEntry(ximp);
    bool flag_veto = false;
    if(fimp->LGFSize()==0 || fimp->LGBSize()==0 || fimp->LGFMax().GetEnergy()<0 || fimp->LGBMax().GetEnergy()<0) continue;
    // both lgf + lgb
    countt++;
    //veto check
    if(fimp->SSSDSize()>0){ flag_veto = veto(fimp->LGFMax(), fimp->SSSD()); }

    // out of the veto
    if(flag_veto) continue;
    count0++;
    pixel piximp = fimp->LGPixel();
    double pin1t = fimp->Pin1T();
    //Clean Map
    CleanMap(&fImpMap, pin1t, 1);
    CleanMap(&fBLMap, pin1t);
    
    //Correlation
    while(xdec<ndec){
      fdec->Clear();
      chdec->GetEntry(xdec++);
      if(fdec->HGFSize()==0 || fdec->HGBSize()==0) continue;
      double dect = fdec->HGFMax().GetTimestamp();
      if(dect>pin1t){
        xdec--;
        break;
      }
      pixel pixdec = DecayPixel(fdec->HGF(), fdec->HGB());
      if(pixdec.first <0 || pixdec.second <0) continue;//don't find any pixel satisify all gates
      pixel corpix = Correlation(fImpMap, pixdec, dect,SA); 
      if(corpix.first<0) continue; // no existing imp correlated to this decay
      ffdec->Clear();
      ffdec->Set(&(fdec->fHits));
      double decaytime = dect - fImpMap[corpix]->Pin1T();
      decaytime = decaytime/1e6; // unit:us
      ffdec->SetDecayTime(decaytime);
      fDecMap[corpix]->push_back(*ffdec);
    }
    

    if(check(fBLMap, piximp, SA)){//reset the bl clock
      Insert(&fBLMap, piximp, fimp->fHits);
      count3++;
    }else{
      std::vector<pixel> vecpix;
      int found = check(fImpMap, piximp, SA, &vecpix);
      if(found){ //back to back
        fImpMap[piximp]->Clear();
        Insert(&fBLMap, piximp, fimp->fHits);
        count2+=1;
        for(int i=0;i<found;i++){
          pixel temp = vecpix[i];
          if(temp.first==piximp.first && temp.second==piximp.second) {count2++;continue;}
          Insert(&fBLMap, temp, fImpMap[temp]->fHits);
          fImpMap[temp]->Clear();
          count2+=1;
        }
      }else{
        Insert(&fImpMap, piximp, fimp->fHits);
      }
    }



    if((ximp%500)==0){
      printf("run%i.root: on entry %lu / %lu \r", runnum,ximp,nimp);
      fflush(stdout);
    }

  }
  
  //clean all fImpMap
  CleanMapAll(&fImpMap);
  printf("run%i.root: on entry %lu / %lu \n", runnum,ximp,nimp);

  std::cout<< "lgf+lgb               = " << countt<< std::endl;
  std::cout<< "out of veto           = " << count0 << std::endl;
  std::cout<< "imp written in TTree  = " << count1 << std::endl;
  std::cout<< "back to back          = " << count2 << std::endl;
  std::cout<< "implanted on blocked  = " << count3 << std::endl;
  std::cout<< "multiple imp correlated= " << count4 << std::endl;

}

