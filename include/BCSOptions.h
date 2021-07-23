
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

  private: 
    bool fSortAndQuit;
    bool fWriteTree;
    bool fWriteEventTree;


  ClassDef(BCSOptions,0)
};

#endif

