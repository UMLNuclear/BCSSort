
#include<iostream>

#include <util.h>

#include <TH1D.h>
#include <BCSint.h>
#include <BCSintFunctions.h>

#include <DDASEvent.h>
#include <BCSOptions.h>

int main(int argc, char **argv){

  BCSint *bcs = BCSint::Get();

  DDASEvent junk1; // trying to force root to load libraries

  for(int i=1;i<argc;i++) {
    //std::cout << "passed " << argv[i] << std::endl;
    std::string sargv = argv[i];
    if(!sargv.substr(sargv.find_last_of(".")+1).compare("root")) { 
      OpenRootFile(argv[i]);
    }
  }

  if(BCSOptions::Get()->SortAndQuit()==true) {
    //bcs->DoSort(); // organize data from run.file to correlated beta.file
    //bcs->TOFfluctuation(); // get TH2D about TOF fluctuation as runtime;
    bcs->CorrectTOF(); // Correct TOF fluctuation. Require tof.root exist;
    bcs->Terminate(0);
  } else {
    bcs->Run(false);
  }


  return 0;
}
