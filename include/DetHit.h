#ifndef DetHit_H
#define DetHit_H
#include <TObject.h>
#include <TH2.h>

class ddaschannel;


class DetHit{
  public: 
    DetHit();
    DetHit(ddaschannel *chan);
    DetHit(const DetHit&);
    ~DetHit();
   

 
    
    void print();/* {
        cout<<number <<"\t"<< energy <<"\t"<< timestamp<<endl;
    }*/
    void SetAddress(int c) {address = c;}   
    void SetNumber(int c) {number = c;}   
    void SetTimestamp(double c) {timestamp = c;}   
    void SetCharge(double c) {charge = c;}   
    int GetAddress()   const { return address; }
    int GetNumber()    const { return number; }
    double GetCharge()    const { return charge; }
    double GetTimestamp() const { return timestamp; }

 
    double GetEnergy() const;
    

  private:    
    int address;
    int number;
    double timestamp;
    double charge;


  ClassDef(DetHit,2)
};


class BCSEvent {
  public:
    BCSEvent();
    virtual ~BCSEvent();

    double Pin1E() const; 
    double Pin1T() const; 
    double Pin2E() const; 
    double Pin2T() const; 
    double Pin3E() const; 
    double Pin3T() const; 
    double I2S() const; //pin1_i2s 
    double PIN1_I2N() const; //pin1_i2n 
    double I2N_I2S() const; 

    std::vector<DetHit> LGF() const;
    std::vector<DetHit> LGB() const;
    std::vector<DetHit> HGF() const;
    std::vector<DetHit> HGB() const;

    int HGFSize() const;
    int HGBSize() const;
    int LGFSize() const;
    int LGBSize() const;
    int SSSDSize() const;

    int  Size() const { return fHits.size(); }
    void Clear() { fHits.clear(); }
    
    int LGPixel() const;
    std::pair<int,int> HGPixel() const;

    TH2D *DrawHG(Option_t *opt="") const;

    bool Range(int x,int low, int high) const { return (x>=low && x<=high); }
    void PrintNumber();

  //private:
 
    std::vector<DetHit> fHits;

  ClassDef(BCSEvent,1)
};

#endif 
