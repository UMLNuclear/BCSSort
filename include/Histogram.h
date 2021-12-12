
#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__

#include <map>
#include <DetHit.h>
#include <Implant.h>

class Histogram {

  private:
    Histogram();
    static Histogram *fHistogram;

  public:
    static Histogram *Get();
    virtual ~Histogram();
   
    std::map<int,double[4][4]> ReadMat(std::string filename = "/home/zhu/notebooks/Ne31/Ne31_Calibration/hpge/xtalmat.dat");
    void EventSort();
    void EventHitPad();
    void ListSort();
    void Process(std::vector<DetHit> vec, std::map<int,double[4][4]> mat);
    void LogX();

  ClassDef(Histogram,0)

};

#endif



















