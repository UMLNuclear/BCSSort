#ifndef __UTIL_H__
#define __UTIL_H__

#include <string>
#include <TCutG.h>
#include <TList.h>


//class TH1;
//class TList;

extern TList *gList;

//void ReadDetMap(std::string filename="detmap.txt");
void ReadDetMap(const char* filename="/home/zhu/packages/BCSSort/config/detmap_10.txt",Option_t *opt="replace");

std::string GetRunNumber(std::string input);

//void ProgressUpdate(long x, long n, bool end = false);

void FillHistogram(std::string hname, TCutG* cut,int xbins, double xlow, double xhigh, double xvalue,
                              int ybins=-1, double ylow=-1, double yhigh=-1, double yvalue=-1);
void FillHistogram(std::string hname,int xbins, double xlow, double xhigh, double xvalue,
                              int ybins=-1, double ylow=-1, double yhigh=-1, double yvalue=-1);
void SaveHistograms(std::string fname="output.root");

#endif
