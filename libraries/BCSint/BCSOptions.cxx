
#include <BCSOptions.h>

BCSOptions *BCSOptions::fOptions = 0;



BCSOptions *BCSOptions::Get() {
  if(!fOptions) 
    fOptions = new BCSOptions;
  return fOptions;
}

BCSOptions::BCSOptions() { Init(); }

BCSOptions::~BCSOptions() { }


void BCSOptions::Init() { //TODO make it read these from a file.
  //fSortAndQuit = false;
  //fWriteTree   = true;
  //fWriteEventTree   = true;
  fSortAndQuit = true;
  fWriteTree   = false;
  fWriteEventTree   = false;

}

