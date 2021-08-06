#ifndef __GH1D_H__
#define __GH1D_H__

#include <TH1D.h>

class GH1D : public TH1D {

  public:
    GH1D();
    GH1D(const TH1D& rhs);
    GH1D(const char *name,int nbins, double low, double high);
    GH1D(const char *name,const char *title, int nbins, double low, double high);
    virtual ~GH1D();

    void Zoom(double low=sqrt(-1),double high=sqrt(-1));

    int Write(const char *name="",Int_t option=0,Int_t bufsize=0);  

  private:

    ClassDef(GH1D,1)
};

#endif
