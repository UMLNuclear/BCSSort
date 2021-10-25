
#ifndef __HISTOGRAM_H__
#define __HISTOGRAM_H__



class Histogram {

  private:
    Histogram();
    static Histogram *fHistogram;

  public:
    static Histogram *Get();
    virtual ~Histogram();
    
    void EventSort();
    void BetaSort();
    void LogX();
    void NaSort();

  ClassDef(Histogram,0)

};

#endif



















