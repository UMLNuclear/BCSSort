
#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>


TGraphErrors* convert_txt_gr(){

  double x[16]    = {0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11, 12 ,13 ,14 ,15}; 
  double xerr[16]; 
  //== pixel 1*1 ==//
  //double y[16] = {3.701,4.212,4.067,4.678,4.923,5.379,4.495,6.021,6.583, 6.283, 5.085, 6.995, 6.506, 6.785, 6.473, 6.579}; 
  //double yerr[16] = {0.484, 0.499, 0.433, 0.43 , 0.417, 0.451, 0.304, 0.475, 0.593, 0.392, 0.259, 0.39 , 0.365, 0, 0.438, 0}; 

  //== pixel 3*3 ==//
  //double y[16] = {4.475,4.332,5.161,5.345,6.385,6.855,5.993,6.928,6.449,6.261,6.25,6.917,7.114,7,5.341,5.437}; 
  //double yerr[16] = {0.453,0,0.393,0.345,0.408,0.441,0.283,0.349,0.296,0.239,0.21,0.222,0.24,0.2,0.17,0.206}; 
 
  double y[16] = {3.921, 5.438, 5.902, 5.621, 7.331, 6.61, 6.542, 6.011, 6.807, 6.36, 6.097, 7.236, 6.271, 6.791, 7.749, 6.936}; 
  double yerr[16] = {0.374, 0.423, 0.395, 0.336, 0.453, 0.43, 0.39, 0.244, 0.327, 0.2, 0.183, 0.236, 0.182, 0.236, 0.293, 0.35}; 
  
  //std::ifstream infile;
  //std::string line;
  //infile.open(filename.c_str());
  //
  //int i=0;
  //while(getline(infile, line)){
  //  if(line[0]=='#') continue;
  //  double cutnum;
  //  double p;
  //  double perr;

  //  std::stringstream ss(line);
  //  ss >> cutnum;
  //  ss >> p;
  //  ss >> perr;

  //  x[i]    = cutnum;
  //  xerr[i] = 0.1;
  //  y[i]    = p;
  //  yerr[i] = perr;
  //  i+=1;

  //}

  TGraphErrors *gr = new TGraphErrors(16,x,y,xerr,yerr);
  return gr;
}
