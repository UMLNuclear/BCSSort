
#ifndef __OUTPUT_MANAGER_H__
#define __OUTPUT_MANAGER_H__

#include <TObject.h>

class TFile;
class TTree;
class Beta;
class Implant;
class Decay;
class BCSEvent;

#include <DetHit.h>

class OutputManager : public TObject {
 
  private:
    OutputManager();
    static OutputManager *fOutputManager;

  public:
    static OutputManager *Get();
    virtual ~OutputManager();


    void Set(std::string runnumber);
    void Fill(Implant *implant,std::vector<Decay> *vdec);
    void FillEvent(std::vector<DetHit> *hits);
    void FillList(DetHit *hit);
    void Close();

    const char *GetName() const;

  private:
    TFile *fFile;
    TTree *fTree;
    Beta  *fBeta;


    TFile *fEventFile;
    TTree *fEventTree;
    BCSEvent  *fEvent;


    TFile *fListFile;
    TTree *fHitTree;
    DetHit  *fHit;

  ClassDef(OutputManager,0)
};

#endif

