
cmake_minimum_required(VERSION 3.10)

project(BCSSort)

find_package(ROOT)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

#set(ARG_NOINSTALL TRUE GLOBAL)
#set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}) 

#set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)


#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib GLOBAL)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


ROOT_GENERATE_DICTIONARY(G__BCSint BCSint.h BCSintFunctions.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/BCSint/LinkDef.h)
add_library(BCSint SHARED G__BCSint.cxx libraries/BCSint/BCSint.cxx libraries/BCSint/BCSintFunctions.cxx)
target_link_libraries(BCSint PUBLIC ${ROOT_LIBRARIES}) 

add_executable(bcssort src/bcssort.cxx) 

target_link_libraries(bcssort ${ROOT_LIBRARIES} BCSint)



