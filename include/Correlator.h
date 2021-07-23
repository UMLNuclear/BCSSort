
#ifndef __CORRELATOR_H__
#define __CORRELATOR_H__


#include <TObject.h>
#include <map>

#include <globals.h>

class Implant;
class Decay;
class DetHit;

class Correlator : public TObject{

  private:
    Correlator();
    static Correlator *fCorrelator;
    Implant *fImplant;
    Decay *fDecay;
    std::map<pixel, Implant*> fImplantMap;
    std::map<pixel, std::vector<Decay> *> fDecayMap;
    std::map<pixel, double> fBlockMap;

    void InitMaps();  

  public:
    static Correlator *Get();
    virtual ~Correlator();

    void AddEvent(std::vector<DetHit> *event);    
    void AddImplant(std::vector<DetHit> *event);
    void AddDecay(std::vector<DetHit> *event);
    pixel Match(pixel dpix, double T);
    void CleanImplantMap(double T);
    void FlushPixel(pixel pix, bool histogram=true);
    void FlushAll(bool histogram=true);

  
  ClassDef(Correlator,0)






};

#endif

