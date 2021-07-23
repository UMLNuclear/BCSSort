
#include <iostream>
#include <string>

int main() {

  std::cout << "hello" << std::endl;

  std::string thing;

  while(true) {
    std::cin >> thing;
    if(thing.compare("quit")==0) 
        break;
    else if(thing.compare("quitfaster")==0) 
        exit(0);
    std::cout << "you wrote: " << thing << std::endl;
  }

  std::cout << "bye, bye" << std::endl;

  return 0;

}




