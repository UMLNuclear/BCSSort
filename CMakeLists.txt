
cmake_minimum_required(VERSION 3.10)

project(BCSSort)

find_package(ROOT)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3 -Wno-unused-variable")
#set(CMAKE_CXX_FLAGS "-g -Wno-unused-variable")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-variable")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

#set(ARG_NOINSTALL TRUE GLOBAL)
#set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}) 

#set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)


#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib GLOBAL)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)                       # .so file. 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


ROOT_GENERATE_DICTIONARY(G__BCSint BCSint.h BCSintFunctions.h BCSOptions.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/BCSint/LinkDef.h)   #create .cxx of root dictionary 
add_library(BCSint SHARED G__BCSint.cxx libraries/BCSint/BCSint.cxx libraries/BCSint/BCSintFunctions.cxx libraries/BCSint/BCSOptions.cxx)               # compile .cxx root dictionary => .o file; "SHARED" = gonna create .so rather than .a file
target_link_libraries(BCSint PUBLIC ${ROOT_LIBRARIES})                                                                  # link all libraries => .so file  


ROOT_GENERATE_DICTIONARY(G__BCSUtil TChannel.h util.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/Util/LinkDef.h)   
add_library(BCSUtil SHARED G__BCSUtil.cxx libraries/Util/TChannel.cxx libraries/Util/util.cxx)               
target_link_libraries(BCSUtil PUBLIC ${ROOT_LIBRARIES})                                                                  


ROOT_GENERATE_DICTIONARY(G__DDAS DDASEvent.h DDASHit.h ddaschannel.h  LINKDEF ${CMAKE_BINARY_DIR}/libraries/DDAS/LinkDef.h)   
#FILE(GLOB DDASSources libraries/DDAS/*.cxx)
add_library(DDAS SHARED G__DDAS.cxx libraries/DDAS/ddaschannel.cxx libraries/DDAS/DDASEvent.cxx libraries/DDAS/DDASHit.cxx libraries/DDAS/DDASHitUnpacker.cxx)               
target_link_libraries(DDAS PUBLIC ${ROOT_LIBRARIES})                                                                  

ROOT_GENERATE_DICTIONARY(G__HIT DetHit.h Implant.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/HIT/LinkDef.h)   
add_library(HIT SHARED G__HIT.cxx libraries/HIT/DetHit.cxx libraries/HIT/Implant.cxx)               

ROOT_GENERATE_DICTIONARY(G__CORRELATOR Correlator.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/Correlator/LinkDef.h)   
add_library(Correlator SHARED G__CORRELATOR.cxx libraries/Correlator/Correlator.cxx)               
target_link_libraries(Correlator PUBLIC ${ROOT_LIBRARIES})                                                                  

ROOT_GENERATE_DICTIONARY(G__OUTPUTMANAGER OutputManager.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/OutputManager/LinkDef.h)   
add_library(OutputManager SHARED G__OUTPUTMANAGER.cxx libraries/OutputManager/OutputManager.cxx)               
target_link_libraries(OutputManager PUBLIC ${ROOT_LIBRARIES})                                                                  

ROOT_GENERATE_DICTIONARY(G__JUNK junk.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/junk/LinkDef.h)   
add_library(Junk SHARED G__JUNK.cxx libraries/junk/junk.cxx)               
target_link_libraries(Junk PUBLIC ${ROOT_LIBRARIES})                                                                  


add_executable(bcssort src/bcssort.cxx) 

target_link_libraries(bcssort ${ROOT_LIBRARIES} BCSint BCSUtil DDAS HIT Correlator OutputManager Junk)



