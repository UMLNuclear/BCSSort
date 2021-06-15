#include <BCSint.h>

#include <pwd.h>

BCSint *BCSint::fInstance=0;

BCSint *BCSint::Get() {
  if(fInstance==0) 
    fInstance = new BCSint();
  return fInstance;
}

BCSint::BCSint(const char *app,int *argc,char **argv,void *options,int numOptions,bool noLogo) 
  : TRint(app,argc,argv,options,numOptions,noLogo) { 

  SetPrompt();

}   

BCSint::~BCSint() { }




const char *BCSint::SetPrompt(const char *newPrompt) {

  return TRint::SetPrompt(newPrompt);
}


void BCSint::Terminate(int status) {
  
  
  //Be polite when you leave.
  printf("\nbye,bye\t%s\n",getpwuid(getuid())->pw_name);
 
  TRint::Terminate(status);
   
}

