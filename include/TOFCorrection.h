
#ifndef __TOFCORRECTION_H__
#define __TOFCORRECTION_H__



class TOFCorrection {

  private:
    TOFCorrection();
    static TOFCorrection *fTOFCorrection;

  public:
    static TOFCorrection *Get();
    virtual ~TOFCorrection();
    
    void Fluctuation();
    void Correct();
    std::vector<double> ReadFile(int num, std::string infilename = "/home/zhu/packages/BCSSort/config/TOF/TOFParameters.txt");

  ClassDef(TOFCorrection,0)

};

#endif



















