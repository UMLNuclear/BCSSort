#ifndef __IMPLANT_H__
#define __IMPLANT_H__

#include<vector>
#include<TObject.h>
#include<TCutG.h>

#include<cstdio>
#include<map>

#include<DetHit.h>

#include <globals.h>

class TH2;
/*
class DetHit{
    public:
        DetHit();
        ~DetHit();


        double GetEnergy() const;

        //private:

        int fNumber;
        unsigned int fAddress;
        double fCharge; //before calibration
        double fTimestamp;

        ClassDef(DetHit,1);
};
*/

/*class HPGeHit:public DetHit{
  public:
  HPGeHit(DetHit);
  ~HPGeHit();


  void Add(HPGeHit);
  }*/




class Implant{
    public:
        Implant();
        Implant(std::vector<DetHit> *hits);
        ~Implant();


        void Clear();
        void Set(std::vector<DetHit> *hits);
        void Print() const; //const(constant fucntion): variables are constant. we cannot change them.
        void PrintFun(TCutG *mycut=0) const;
        void SimplePrint() const;
        TH2 *Draw(); // return TH2D

        int FrontSize() const { return fDSSDFront.size(); }
        int BackSize() const { return fDSSDBack.size(); }
        bool Stopped() const; // {return !fSSSD.size(); } 
        bool Inside(TCutG *mycut) {return mycut->IsInside(fI2S, fPIN1E);} 

        bool IsGood() const ;
        std::pair<int,int> GetPixel() const;   
        void Copy(Implant &) const;     

        double GetTimestamp() const { return fPIN1T; }

        //private:
        double fPIN1E;
        double fPIN1T;

        double fPIN2E;
        double fPIN2T;

        double fI2N;
        double fI2S;

        double fI2NT;
        double fI2ST;

        double fI2S_I2N;
        double fI2S_I2N_T;

        std::vector<DetHit> fDSSDFront;
        std::vector<DetHit> fDSSDBack;

        std::vector<DetHit> fSSSD;

        ClassDef(Implant,1);
};


class Decay{
    public:
        Decay();
        ~Decay();

        void Set(std::vector<DetHit> *hits);
        void Clear();

        std::vector<DetHit> fDSSDFront;
        std::vector<DetHit> fDSSDBack;
        std::vector<DetHit> fLaBr;
        std::vector<DetHit> fGe;
        std::vector<DetHit> fSSSD;
        double fImplantTime;

        int FrontSize() const { return fDSSDFront.size(); }
        int BackSize() const { return fDSSDBack.size(); }
        int GeSize() const { return fGe.size(); }
        void SimplePrint(double) const;

        void SetImplantTime(double impt) {fImplantTime = impt;}//pin1.t - fh.t
        //bool Stopped() const {return !fSSSD.size(); } 

        double GetTimestamp() const {
          if(FrontSize()<1) return -1;
          return fDSSDFront.front().GetTimestamp();
        }

        
        bool IsGood() const;        
        std::pair<int,int> GetPixel() const;   
        

        ClassDef(Decay,1);
};





class Clover{
    public:
        Clover(){}
        ~Clover(){}


        std::vector<DetHit> fXtal;
        //int fID;
        double fAddE; 


        void Clear(){
            fXtal.clear();
            //fID = -1;
        }
        void Add(DetHit hit) {fXtal.push_back(hit);}

        int Size() const {return fXtal.size();}
        double AddbackSum() const {   //// addback energy of each clover in each event
            double sum = 0;
            for (size_t x=0;x<Size();x++){
                sum += fXtal[x].GetEnergy();
            }
            return sum;
        }



        std::vector<double> AddbackSum(TCutG *gate) const;   // for summing

        void SetAddE(double e=-1){
            if(e<0){
                fAddE=AddbackSum();
            }else{
                fAddE = e;
            }
        }

        double AddbackE() const {return fAddE;} // must be called after "SetAddE()"
        ClassDef(Clover,1);

};

class Beta{
  public:
    Beta(){}
    Beta(Implant &imp, std::vector<Decay> &dec){fImplant = imp; fDecay = dec;}
    ~Beta(){}

    
    void Set(Implant &imp, std::vector<Decay> &dec){fImplant = imp; fDecay = dec;}

    int DecaySize() const { return fDecay.size(); }
    void Clear(); 
    void SimplePrint() const;
    bool Inside(TCutG *mycut) {return fImplant.Inside(mycut);}   
    
    bool Stopped() const {return !fImplant.fSSSD.size(); } 

    pixel GetImplantPixel() const { return fImplant.GetPixel(); }


  //private: 
    Implant fImplant;
    std::vector<Decay> fDecay;

    ClassDef(Beta,1);

};



#endif
