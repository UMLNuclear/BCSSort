#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <TCutG.h>
#include <TList.h>


//class TH1;
//class TList;

extern TList *gList;


std::string GetRunNumber(std::string input);


void FillHistogram(std::string hname, TCutG* cut1, TCutG *cut2, int xbins, double xlow, double xhigh, double xvalue,
                              int ybins=-1, double ylow=-1, double yhigh=-1, double yvalue=-1, double zvalue=sqrt(-1));
void FillHistogram(std::string hname, TCutG* cut,int xbins, double xlow, double xhigh, double xvalue,
                              int ybins=-1, double ylow=-1, double yhigh=-1, double yvalue=-1,  double zvalue=sqrt(-1));
void FillHistogram(std::string hname,int xbins, double xlow, double xhigh, double xvalue,
                              int ybins=-1, double ylow=-1, double yhigh=-1, double yvalue=-1, double zvalue=sqrt(-1));
void SaveHistograms(std::string fname="output.root", Option_t *opt="");





#endif
