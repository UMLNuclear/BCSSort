
cmake_minimum_required(VERSION 3.10)

project(BCSSort)

find_package(ROOT)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

file(MAKE_DIRECTORY bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


add_executable(bcssort src/bcssort.cxx) 

target_link_libraries(bcssort ${ROOT_LIBRARIES})
