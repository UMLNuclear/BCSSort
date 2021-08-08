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

    void DoSort(); // organize data
    void TOFfluctuation(); // check tof fluctuation of each event file.
    void CorrectTOF();

    int UpdateFileCount() { fFileCount++; return fFileCount-1; }


  ClassDef(BCSint,0)

};

#endif
