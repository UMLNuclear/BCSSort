
#include <junk.h>


tom *tom::fInstance = 0;



tom* tom::Get() {
    if(fInstance==0) {
      fInstance= new tom();
    }
    return fInstance;
}

