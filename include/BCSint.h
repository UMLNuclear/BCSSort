#ifndef __BCSINT_H__
#define __BCSINT_H__

#include<string>
#include<iostream>

#include <TRint.h>
#include <TChain.h>
#include <TList.h>

extern TChain *gChain;
extern TList  *gCuts;

class BCSint : public TRint { 
  public:
    static BCSint *Get();
    virtual ~BCSint();

  private:
    static BCSint *fInstance;
    BCSint(const char *app="bcssort",int *argc=0,char **argv=0,void *options=0,int numOptions=0,bool noLogo=true);
    int fFileCount;

  public:
    virtual const char *SetPrompt(const char *newPrompt="bcs [%d] ");
    virtual void        Terminate(int status);

    //void OpenRootFile(std::string fname); 

    void DoSort(); // organize data: dchan->beta & event;
    void EventSort(); // organize data: dchan->DetHit
    void Correlation(); //organize data: input(event)-> call(implant&decay)-> build beta;
    void BetaSort(); // from beta to beta: by beta_tof correction.
    void QuickBetaSort();

    int UpdateFileCount() { fFileCount++; return fFileCount-1; }


  ClassDef(BCSint,0)

};

#endif
