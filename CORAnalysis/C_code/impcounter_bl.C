
#define SA 1


typedef std::pair<int, int> pixel;

TCutG *ban;
TCutG *cutstrip;

int countt = 0;
int count0 = 0; // counts out of veto;
int count1 = 0; // good implant will be written in TTree;
int count2 = 0; // black list: back to back;
int count3 = 0; // black list: pixel is still cooling down;

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
void InitMap(std::map<pixel, BCSEvent*> *mp1, std::map<pixel, BCSEvent*> *mp2){
  std::map<pixel, BCSEvent*> temp1;
  std::map<pixel, BCSEvent*> temp2;
  for(int m=0;m<40;m++){
    for(int n=0;n<40;n++){
      pixel pix = std::make_pair(m,n);
      temp1[pix] = new BCSEvent;
      temp2[pix] = new BCSEvent;
    }
  }
  *mp1 = temp1;
  *mp2 = temp2;
}

//=========== Clean Map ===============//
void CleanMap(std::map<pixel, BCSEvent*> *mp, double pin1t, int flush=0){
  for(auto &it: (*mp)){
    if(it.second->Pin1T()<0) continue;
    double dt = pin1t - it.second->Pin1T();
    if(dt>1e9){
      if(flush) count1++;
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
void insert(std::map<pixel, BCSEvent*> *mp, pixel pix, std::vector<DetHit> hits){
  (*mp)[pix]->Clear();
  (*mp)[pix]->fHits = hits;
}

//================== "main" code for running ===============//
void counter(int runnum=1040){

  TChannel::ReadDetMapFile();
  TFile *cutf2 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/bananagate2.root");
  ban = (TCutG *)cutf2->Get("ban");

  TFile *cutf3 = TFile::Open("/home/zhu/packages/BCSSort/root_file/new_cut/stripgate.root");
  cutstrip = (TCutG *)cutf3->Get("cutstrip");

  TChain *chan = new TChain("event");
  chan->Add(Form("/home/zhu/packages/BCSSort/CORAnalysis/example/data/implant/event%i*root",runnum));
  BCSEvent *fevent = new BCSEvent;
  chan->SetBranchAddress("BCSEvent", &fevent);
  long nentries = chan->GetEntries();
  long x = 0;

  std::map<pixel, BCSEvent*>fImpMap;
  std::map<pixel, BCSEvent*>fBLMap;

  //Initialize map
  InitMap(&fImpMap, &fBLMap); 

  for(x=0;x<nentries;x++){
    fevent->Clear();
    chan->GetEntry(x);
    bool flag_veto = false;
    if(fevent->LGFSize()==0 || fevent->LGBSize()==0 || fevent->LGFMax().GetEnergy()<0 || fevent->LGBMax().GetEnergy()<0) continue;
    // both lgf + lgb
    countt++;
    //veto check
    if(fevent->SSSDSize()>0){ flag_veto = veto(fevent->LGFMax(), fevent->SSSD()); }

    // out of the veto
    if(flag_veto) continue;
    count0++;
    pixel piximp = fevent->LGPixel();
    double pin1t = fevent->Pin1T();
    //Clean Map
    CleanMap(&fImpMap, pin1t,1);
    CleanMap(&fBLMap, pin1t);
    //Map Cleaned
    if(check(fBLMap, piximp, SA)){//reset the bl clock
      insert(&fBLMap, piximp, fevent->fHits);
      count3++;
    }else{
      std::vector<pixel> vecpix;
      int found = check(fImpMap, piximp, SA, &vecpix);
      if(found){ //back to back
        fImpMap[piximp]->Clear();
        insert(&fBLMap, piximp, fevent->fHits);
        count2+=1;
        for(int i=0;i<found;i++){
          pixel temp = vecpix[i];
          if(temp.first==piximp.first && temp.second==piximp.second) {count2++;continue;}
          insert(&fBLMap, temp, fImpMap[temp]->fHits);
          fImpMap[temp]->Clear();
          count2+=1;
        }
      }else{
        insert(&fImpMap, piximp, fevent->fHits);
      }
    }



    if((x%5000)==0){
      printf("run%i.root: on entry %lu / %lu \r", runnum,x,nentries);
      fflush(stdout);
    }

  }
  
  //clean all fImpMap
  CleanMapAll(&fImpMap);
  printf("run%i.root: on entry %lu / %lu \n", runnum,x,nentries);

  std::cout<< "lgf+lgb               = " << countt<< std::endl;
  std::cout<< "out of veto           = " << count0 << std::endl;
  std::cout<< "imp written in TTree  = " << count1 << std::endl;
  std::cout<< "back to back          = " << count2 << std::endl;
  std::cout<< "implanted on blocked  = " << count3 << std::endl;

}

