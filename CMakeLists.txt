
cmake_minimum_required(VERSION 3.10)

project(BCSSort)

file(MAKE_DIRECTORY bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(bcssort src/bcssort.cxx) 


