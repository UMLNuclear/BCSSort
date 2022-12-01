#include<util.h>
#include<globals.h>

void implantloop(){

  TChannel::ReadDetMapFile();


  std::map<pixel, double>  fImpTMap;
  std::map<pixel, Implant*>fImpMap;

  int dnum;
  double dt, dt_nb;



  int c0 = 0;
  int c1 = 0;
  int c2 = 0;
  int c3 = 0;
  int c4 = 0;

  
  int runnum = 1032;
  TChain *imp = new TChain("implant");
  imp->Add(Form("/home/zhu/packages/BCSSort/data/5us_tofcor/implant/implant%i*.root",runnum));
  long entries = imp->GetEntries();
  Implant *fimp = new Implant;
  imp->SetBranchAddress("Implant",&fimp);
  
  
  long x = 0;
  for(x=0;x<entries;x++){
    fimp->Clear();
    imp->GetEntry(x);
    if(fimp->FrontSize()>2) c1++;    
    if(fimp->BackSize()>2) c2++;    
 
   
    if((x%5000)==0){
      printf("on entry %lu / %lu \r",x,entries);
      fflush(stdout);
    } 
  }


  printf("on entry %lu / %lu \n\n\n",x,entries);


  std::cout << "total entries           = " << entries << std::endl;
  std::cout << "LGF Fired over 2 strips = " << c1      << std::endl;
  std::cout << "LGB Fired over 2 strips = " << c2      << std::endl;


}
