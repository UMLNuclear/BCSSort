


#include <util.h>

#include <TH1D.h>

#include <BCSint.h>

int main(int argc, char **argv){

  //Hello();
  //printf("yo\n");
  //TH1D h("h","h",100,0,100);
  BCSint *bcs = BCSint::Get();
  bcs->Run(1);
  
  return 0;
}
