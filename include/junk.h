
#ifndef __JUNK_H__
#define __JUNK_H__


#include<string>


class cat {
  public:
      cat() {}
      ~cat() {}

      std::string sound = "meow";
};

class tom : public cat {
  private:  
    tom() : cat()  { }
    static tom *fInstance;

  public:
      ~tom() { }

      static tom* Get();

      std::string name ="tom";


};


#endif
