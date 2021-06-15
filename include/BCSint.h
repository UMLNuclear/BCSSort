#ifndef __BCSINT_H__
#define __BCSINT_H__

#include <TRint.h>

class BCSint : public TRint { 
  public:
    static BCSint *Get();
    virtual ~BCSint();

  private:
    static BCSint *fInstance;
    BCSint(const char *app="bcssort",int *argc=0,char **argv=0,void *options=0,int numOptions=0,bool noLogo=true);


  public:
    virtual const char *SetPrompt(const char *newPrompt="bcs [%d] ");
    virtual void        Terminate(int status);


  ClassDef(BCSint,0)

};

#endif
