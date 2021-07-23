
#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include <utility>

#define BUILDTIME 5000  //time  in nanoseconds
#define EXPIRETIME 500E6  //time  in nanoseconds

typedef std::pair<int,int> pixel;



#define RESET_COLOR "\033[m"
#define BLUE       "\033[1;34m"
#define YELLOW     "\033[1;33m"
#define GREEN      "\033[1;32m"
#define RED        "\033[1;31m"
#define BLACK      "\033[1;30m"
#define MAGENTA    "\033[1;35m"
#define CYAN       "\033[1;36m"
#define WHITE      "\033[1;37m"

#define DBLUE      "\033[0;34m"
#define DYELLOW    "\033[0;33m"
#define DGREEN     "\033[0;32m"
#define DRED       "\033[0;31m"
#define DBLACK     "\033[0;30m"
#define DMAGENTA   "\033[0;35m"
#define DCYAN      "\033[0;36m"
#define DWHITE     "\033[0;37m"

#define BG_WHITE   "\033[47m"
#define BG_RED     "\033[41m"
#define BG_GREEN   "\033[42m"
#define BG_YELLOW  "\033[43m"
#define BG_BLUE    "\033[44m"
#define BG_MAGENTA "\033[45m"
#define BG_CYAN    "\033[46m"




#endif
