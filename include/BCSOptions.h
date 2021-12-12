
#ifndef __BCSOPTIONS_H__
#define __BCSOPTIONS_H__

#include<TObject.h>

class BCSOptions : public TObject {

  private:
    BCSOptions();
    static BCSOptions *fOptions;
    void Init();


  public:
    static BCSOptions *Get();
    virtual ~BCSOptions();


    bool SortAndQuit() { return fSortAndQuit; }
    bool WriteTree()   { return fWriteTree;   }
    bool WriteEventTree()   { return fWriteEventTree;   }
    bool WriteListTree()   { return fWriteListTree;   }
    bool WriteImpTree()   { return fWriteImpTree;   }
    bool WriteDecTree()   { return fWriteDecTree;   }

  private: 
    bool fSortAndQuit;
    bool fWriteTree;
    bool fWriteEventTree;
    bool fWriteListTree;
    bool fWriteImpTree;
    bool fWriteDecTree;


  ClassDef(BCSOptions,0)
};

#endif

