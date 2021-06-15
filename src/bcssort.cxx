


#include <util.h>

#include <TH1D.h>

#include <BCSint.h>

int main(int argc, char **argv){

  BCSint *bcs = BCSint::Get();
  bcs->Run(false);
 
  return 0;
}
