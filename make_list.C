#include <Implant.h>
#include <DetHit.h>
#include <map>
#include <util.h>
#include <ddaschannel.h>
#include <DDASEvent.h>
#include <TChannel.h>
#include <BCSint.h>
#include <TCutG.h>
#include <TFile.h>



Implant imp;
Decay de;
std::map<Implant*,std::vector<Decay>* > fCoorMap;
std::map<pair<int,int>, Implant*> fPixelMap;
std::map<pair<int,int>, double> fBlockMap;

long decays_with_no_implant = 0;
Beta beta;


vector<TCutG*> cuts;

#define BUILDTIME 5000  //time  in nanoseconds
#define EXPIRETIME 500E6  //time  in nanoseconds


//================ Load TCutG ===================//
void LoadCuts(){
    TFile *mycut = TFile::Open("~zhu/notebooks/Ne31/Time_Correlation/mynewcuts.root");
    TIter keys (mycut->GetListOfKeys());
    while(TKey *key = (TKey*)keys.Next()) {
        cuts.push_back((TCutG*)key->ReadObj());
    }
}

//===================== Beta ===========================//
void HandleBeta(Beta *beta){
    if(!beta) return ;
    FillHistogram("pid_beta", 2e3,0,2e4,beta->fImplant.fI2S, 8e3,0,8e3,beta->fImplant.fPIN1E);
    if(beta->fImplant.Stopped()){
        FillHistogram("pidnosssd_beta", 2e3,0,2e4,beta->fImplant.fI2S, 8e3,0,8e3,beta->fImplant.fPIN1E);
    }
    delete beta;
}


//=================== BlackList =====================//  it could be written as one "if()" in "HandelVecHit()"
bool BlackList(Implant cur, Implant temp){
    if(fPixelMap.count(temp.GetPixel())){
        double dt = fabs(cur.fPIN1T-temp.fPIN1T);
        if(dt<EXPIRETIME){
            return false;   // if cur.pixel = imp.pixel && dt < EXPIRETIME => false => fill this pixel and fPIN1T in fBlockMap
        }else return true;
    }else return true;
}

bool CheckBlockMap(Implant temp){   // could be written as two "if()" in "HandelVecHit()"
    std::pair<int,int> tempixel = temp.GetPixel();
    if(fBlockMap.count(tempixel)){
        double t = fBlockMap[tempixel];
        if(fabs(t - temp.fPIN1T)<EXPIRETIME){
            fBlockMap[tempixel] = imp.fPIN1T;
            return true;            // if temp.pixel is counted in fBlockMap && dt<EXPIRETIME => true => reset time
        }else{
            fBlockMap.erase(tempixel); // if dt>EXPIRETIME || it isn't counted in fBlockMap => false => erase pixel from fBlockMap;
            return false;
        }        
    }else return false;    
}

//===================== Clean Implant Map ====================//
void CleanMap(double t){
    std::map<pair<int,int>, Implant*>::iterator it;
    for(it=fPixelMap.begin();it!=fPixelMap.end();it++){
        double dt = fabs(t - it->second->fPIN1T);
        if(dt>EXPIRETIME){
            Implant *cur = it->second;
            FillHistogram("pid_cur", 2e3,0,2e4,cur->fI2S, 8e3,0,8e3,cur->fPIN1E);
            if(cur->Stopped()){
                FillHistogram("pidnosssd_cur", 2e3,0,2e4,cur->fI2S, 8e3,0,8e3,cur->fPIN1E);
            }
            // TODO: save implant
            //beta.Clear();
            //beta.fImplant = *cur;
            //beta.fDecay = *fCoorMap[cur];
            cur->Clear();
            fCoorMap[cur]->clear();
        }
    }
}

//==================== Search Implant Map =====================//
std::pair<int,int> SearchImpMap(std::pair<int,int> pixel, double t){
    if(fPixelMap[pixel]->IsGood()){
        // TODO: blacklist pixel for sometime(eg, EXPIRETIME), implant for Expiretime with no decay?
        return pixel;
    }else{
        double dt = 0;
        int i = -1;
        int j = -1;
        for(int row=pixel.first-1;row<=pixel.first+1;row++){
            for(int col=pixel.second-1;col<=pixel.second+1;col++){
                if((row<0||row>39) || (col<0||col>39)) continue;
                Implant *neighbor = fPixelMap[std::make_pair(row,col)];
                if(neighbor->IsGood()){
                    if(fabs(t - neighbor->fPIN1T)>dt){
                        dt = fabs(t - neighbor->fPIN1T);
                        i = row;
                        j = col;
                    }
                }
            }
        }
        if(i>-1 && j>-1) return std::make_pair(i,j);
        else return pixel;
    }
}

//=================== Search Correlated Decay ==================//
std::pair<int,int> SearchDecMap(std::pair<int,int> pixel, double t){
    double dt = DBL_MAX;
    std::pair<int,int> depixel = std::make_pair(-1,-1);
    for(int row=pixel.first-1;row<=pixel.first+1;row++){
        for(int col=pixel.second-1;col<=pixel.second+1;col++){
            if((row<0||row>39) || (col<0||col>39)) continue;
            std::pair<int,int> temp = std::make_pair(row,col);
            //if(fPixelMap[temp]->IsGood()){
                if(fabs(t - fPixelMap[temp]->fPIN1T) < dt){
                    dt = t - fPixelMap[temp]->fPIN1T; 
                    depixel = temp;
                }                    
            //}
        }
    }
    return depixel;
}



//==================== Handel Hit Vector ==================//
void HandleVecHit(std::vector<DetHit> *vec){
    bool pin1 = false;
    for(int x=0;x<vec->size();x++){
        if(vec->at(x).GetNumber() == 181){
            pin1 = true;
            break;
        }
    }  
    // Implant //
    if(pin1){
        imp.Set(vec);
        if(imp.IsGood()){
            std::pair<int,int> pixel = imp.GetPixel();
            CleanMap(imp.fPIN1T);
            if(!CheckBlockMap(imp)) {
                Implant *cur = fPixelMap[SearchImpMap(pixel, imp.fPIN1T)];
                if(cur->IsGood()){
                    // BlackList
                    if(BlackList(*cur,imp)){
                        // TODO: fill beta
                        Beta *beta = new Beta(*cur, *fCoorMap[cur]);
                        HandleBeta(beta);
                        //beta.Clear();
                        //beta.fImplant = *cur;
                        //beta.fDecay = *fCoorMap[cur];
                        FillHistogram("pid_cur", 2e3,0,2e4,cur->fI2S, 8e3,0,8e3,cur->fPIN1E);
                        if(cur->Stopped()){
                            FillHistogram("pidnosssd_cur", 2e3,0,2e4,cur->fI2S, 8e3,0,8e3,cur->fPIN1E);
                        }

                        //    // TODO: Decay has no copy constructor, need to fix
                        //    Decay temp = fCoorMap[cur]->at(x);
                        //    std::cout << "\t Decay(" 
                        //              << temp.GetPixel().first
                        //              << ":"<<temp.GetPixel().second 
                        //              << ") @ "
                        //              << (temp.fDSSDFront[0].GetTimestamp() - cur->fPIN1T)/1e6
                        //              << " ms" <<std::endl;
                        //std::cout<< "======================================" <<std::endl;
                        FillHistogram("decay_size",50,0,50,fCoorMap[cur]->size());
                        cur->Clear();
                        fCoorMap[cur]->clear();
                        imp.Copy(*cur);
                    }else{
                        fBlockMap[pixel] = imp.fPIN1T;
                        cur->Clear();
                        fCoorMap[cur]->clear();
                    }
                }else{
                    imp.Copy(*cur);  // copy the imp to cur(should be empty)
                }
            }
        } else {  // BAD IMPLANT EVENTS
            FillHistogram("badimplant_size",40,0,40,imp.FrontSize(),
                    40,0,40,imp.BackSize()); 
        }
    }else{  // Decay events
        de.Set(vec);
        if(de.IsGood()){
            std::pair<int,int> pixel = de.GetPixel();
            double fhe = de.fDSSDFront[0].GetEnergy();
            double bhe = de.fDSSDBack[0].GetEnergy();
            std::pair<int,int> depixel = SearchDecMap(pixel, de.fDSSDFront[0].GetTimestamp());
            if(depixel.first<0 || depixel.second<0){
                decays_with_no_implant++;
            }else if(fabs(fhe-bhe)<((fhe+bhe)/2)*0.05){
                Implant *cur = 0;
                cur = fPixelMap[depixel];
                fCoorMap[cur]->push_back(de);
                double FHtime = de.fDSSDFront[0].GetTimestamp();
                map<int,Clover> clmap;
                for(int y=0;y<de.GeSize();y++){
                    double deltatime = FHtime - de.fGe[y].GetTimestamp();
                    FillHistogram("tdif_gamma_decay", 200,-1000,1000, deltatime, 300,0,300,de.fGe[y].GetNumber());
                    if(deltatime>40 && deltatime<300){ // Fill gamma spectrum
                        if(de.fGe[y].GetEnergy()>10 && de.fGe[y].GetEnergy()<4000){
                            int clnum = (de.fGe[y].GetNumber()-208)/4;
                            clmap[clnum].Add(de.fGe[y]);
                            FillHistogram("summary", 300,0,300,de.fGe[y].GetNumber(), 8e3,0,4e3,de.fGe[y].GetEnergy());
                        }
                    }
                }
                for(int i=0;i<cuts.size();i++){
                    if(cuts[i]->IsInside(cur->fI2S,cur->fPIN1E)){
                        map<int,Clover>::iterator it;
                        for(it=clmap.begin();it!=clmap.end();it++){
                            it->second.SetAddE();
                            FillHistogram(Form("addback_%s",cuts[i]->GetName()), 8e3,0,4e3,it->second.AddbackE());
                            for(size_t m=0;m<it->second.Size();m++){
                                FillHistogram(Form("single_%s", cuts[i]->GetName()), 8e3,0,4e3,it->second.fXtal[m].GetEnergy());
                                for(size_t n=m+1;n<it->second.Size();n++){
                                    FillHistogram(Form("ggmat_%s", cuts[i]->GetName()), 8e3,0,4e3,it->second.fXtal[m].GetEnergy(), 
                                            8e3,0,4e3, it->second.fXtal[n].GetEnergy());
                                    FillHistogram(Form("ggmat_%s", cuts[i]->GetName()), 8e3,0,4e3,it->second.fXtal[n].GetEnergy(), 
                                            8e3,0,4e3, it->second.fXtal[m].GetEnergy());
                                }     
                            }
                        }
                    }
                }
            }
        }
    }
    imp.Clear();
    de.Clear(); 

}

bool PixelCheck(std::pair<int,int> x, std::pair<int,int> y){

    if(((y.first>=(x.first-1)) && (y.first<=(x.first+1))) && 
            ((y.second>=(x.second-1)) && (y.second<=(x.second+1)))  )  {

        return true; 

    }
    return false;
}

void InitCoorMap() {
    for(int x=0;x<40;x++){
        for(int y=0;y<40;y++){
            std::pair<int,int> pixel = std::make_pair(x,y);
            fPixelMap[pixel] = new Implant;
            fCoorMap[fPixelMap[pixel]] = new std::vector<Decay>;    
        }
    }    
}

void CloseCoorMap() {
    for(int x=0;x<40;x++){
        for(int y=0;y<40;y++){
            std::pair<int,int> pixel = std::make_pair(x,y);
            delete fCoorMap[fPixelMap[pixel]];    
            delete fPixelMap[pixel];
        }
    }    
}

void make_list(){
    InitCoorMap();
    LoadCuts();
    DDASEvent *event  = new DDASEvent;
    ddaschannel *chan = new ddaschannel;

    gChain->SetBranchAddress("ddasevent",&event);

    long nentries = gChain->GetEntries();
    //nentries = 1e6;
    long x = 0;

    double starttime = 0;
    double buildtime = BUILDTIME;
    std::vector<DetHit> *vec_hit = new std::vector<DetHit>;

    Beta beta;
    TTree *beta_tree = new TTree("beta_tree","beta_tree");
    beta_tree->Branch("beta", &beta);
    map<pair<int,int>, Beta> pixel_map;
    TChannel::ReadDetMapFile();
    for(x=0;x<nentries;x++) {
        gChain->GetEntry(x);
        //FillHistogram("NEvents",100,0,100,event->GetNEvents());
        for(int y=0;y<event->GetNEvents();y++) {
            chan = event->GetData()[y];          
            DetHit hit(chan);
            //FillHistogram("summary",16000,0,16000,hit.energy,
            //        300,0,300,hit.GetNumber());
            //chan->Print();



            if((hit.GetTimestamp()-starttime)>buildtime){
                HandleVecHit(vec_hit);
                starttime = hit.GetTimestamp();
                vec_hit->clear();
            }
            vec_hit->push_back(hit);


            /**************************** Event ******************************/
            /*if((hit.GetTimestamp()-starttime)>buildtime){
              starttime = hit.GetTimestamp();
              if(hits->size()>0){
              imp.Set(hits);
              if(imp.fPIN1T<=0){
              de.Set(hits);
              imp.Clear();
              }
              }
              hits->clear();
              int indexf = -1;
              int indexb = -1;
              if(imp.FrontSize() && imp.BackSize()){
              indexf = 0;
              indexb = 0;
              if(imp.FrontSize()>1){
              double tempE = imp.fDSSDFront[0].GetEnergy();
              for(int z=1;z<imp.FrontSize();z++){
              if(tempE<imp.fDSSDFront[z].GetEnergy()){
              tempE = imp.fDSSDFront[z].GetEnergy();
              indexf = z;
              }
              }
              }
              if(imp.BackSize()>1){
              double tempE = imp.fDSSDBack[0].GetEnergy();
              for(int z=1;z<imp.BackSize();z++){
              if(tempE<imp.fDSSDBack[z].GetEnergy()){
              tempE = imp.fDSSDBack[z].GetEnergy();
              indexb = z;
              }
              }
              }
              }
              if(indexf>=0 && indexb>=0){
              int f = imp.fDSSDFront[indexf].GetNumber();
              int b = imp.fDSSDBack[indexb].GetNumber();
              pair<int,int> pixel = make_pair(f,b); // new implant
              map<pair<int,int>, Beta>::iterator it;
              for(it=pixel_map.begin();it!=pixel_map.end();it++){
              double dt = fabs(imp.fPIN1T - it->second.fImplant.fPIN1T);
              if(dt>500e6){
              beta.Clear();
              beta = it->second;
              beta_tree->Fill();
              std::cout<<beta_tree->GetEntries()<<std::endl;
              pixel_map[it->first].Clear();
              pixel_map.erase(it->first);
              }
              }
              bool found = 0;
              int i = -1;
              int j = -1;
            //pair<int,int> replace_pixel;
            if(pixel_map.count(pixel)){
            found = 1;
            pair<int,int> replace_pixel = pixel;
            }else{
            double dt_pixel = 0;
            for(int row=pixel.first-1;row<=pixel.first+1;row++){
            for(int col=pixel.second-1;col<=pixel.second+1;col++){
            pair<int,int>temp_pixel = make_pair(row,col);
            if(pixel_map.count(temp_pixel)){
            if(fabs(imp.fPIN1T-pixel_map[temp_pixel].fImplant.fPIN1T)>dt_pixel){
            dt_pixel = fabs(imp.fPIN1T-pixel_map[temp_pixel].fImplant.fPIN1T);
            i = row;
            j = col;
            found = 1; 
            }    
            }
            }
        }
        if(i>=0 && j>=0){
            pair<int,int> replace_pixel = make_pair(i,j);
        }
        }
        if(!found){
            Beta temp;
            temp.fImplant = imp;
            pixel_map[pixel] = temp;
        }
        if(found){
            beta.Clear();
            beta = pixel_map[replace_pixel];
            beta_tree->Fill();
            pixel_map[replace_pixel].Clear();
            pixel_map.erase(replace_pixel);
            Beta temp;
            temp.fImplant = imp;
            pixel_map[pixel] = temp;
        }
        }
        if(de.FrontSize()==1 && de.BackSize()==1){
            int i = -1;
            int j = -1;
            int f = de.fDSSDFront[0].GetNumber() + 40;
            double time = de.fDSSDFront[0].GetTimestamp();
            double fhe = de.fDSSDFront[0].GetEnergy();
            double bhe = de.fDSSDBack[0].GetEnergy();
            int b = de.fDSSDBack[0].GetNumber() + 40;
            double dtime = DBL_MAX;
            if(fabs(fhe-bhe)<((fhe+bhe)/2)*0.05){
                for(int row=f-1;row<=f+1;row++){
                    for(int col=b-1;col<=b+1;col++){
                        pair<int,int> pixel_temp = make_pair(row,col);
                        if(pixel_map.count(pixel_temp)){
                            if(fabs(time-pixel_map[pixel_temp].fImplant.fPIN1T)<dtime){
                                dtime = time - pixel_map[pixel].fImplant.fPIN1T;
                                i = row;
                                j = col;
                            }
                        }
                    }
                }
            }
            if(i<=0 || j<=0){

            }else{
                pair<int,int> pixel = make_pair(i,j);
                //de.SetImplantTime(de.fDSSDFront[0].GetTimestamp()-pixel_map[pixel].fImplant.fPIN1T);
                pixel_map[pixel].fDecay.push_back(de);
            }
        }

        }
        hits->push_back(hit);*/
        }
        if((x%50000)==0){
            printf("  on entry %lu / %lu            \r", x,nentries);
            fflush(stdout);
        } 
    }

    printf("  on entry %lu / %lu            \n", x,nentries);
    CloseCoorMap();
}
