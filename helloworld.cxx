


#include<iostream>
#include<string>


enum animal{
  cat = 100, 
  dog = 1000,
  mouse = 7
};

std::string speak(animal a) {
  std::string sound=""; 
  switch(a) {
    case cat:
      sound = "meow";
      break;
    case dog:
      sound = "bark";
      break;
    case mouse:
      sound = "squek";
      break;
  }
  return sound;
}



int main(){



  std::cout<<"hello"<<std::endl;
  std::cout<<"the cat goes " << speak(cat) <<std::endl;
  std::cout<<"the dogt goes " << speak(dog) <<std::endl;
  std::cout<<"the mouse goes " << speak(mouse) <<std::endl;



  return 0;
}
