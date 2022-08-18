
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
   
    std::map<int,double[4][4]> ReadMat(std::string filename = "/home/zhu/packages/BCSSort/config/ct_correction.dat");
    void EventSort();
    void BetaSort();
    void Beta150GateTOF();
    void Beta3DPID();
    void BetaPID();
    void BetaTOF();
    void ListSort();
    void Process(std::vector<DetHit> vec, std::map<int,double[4][4]> mat);
    void LogX();

  ClassDef(Histogram,0)

};

#endif



















