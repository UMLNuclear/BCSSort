
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
  fSortAndQuit = false;
  fWriteTree   = false;
  fWriteEventTree   = false;
  fWriteListTree   = false;
  fWriteImpTree   = false;
  fWriteDecTree   = false;
  
  fSortAndQuit = true;
  //fWriteTree   = true;
  fWriteEventTree   = true;
  fWriteListTree   = true;
  //fWriteImpTree   = true;
  //fWriteDecTree   = true;

}

