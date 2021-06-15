
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
set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)                       # .so file. 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


ROOT_GENERATE_DICTIONARY(G__BCSint BCSint.h BCSintFunctions.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/BCSint/LinkDef.h)   #create .cxx of root dictionary 
add_library(BCSint SHARED G__BCSint.cxx libraries/BCSint/BCSint.cxx libraries/BCSint/BCSintFunctions.cxx)               # compile .cxx root dictionary => .o file; "SHARED" = gonna create .so rather than .a file
target_link_libraries(BCSint PUBLIC ${ROOT_LIBRARIES})                                                                  # link all libraries => .so file  


ROOT_GENERATE_DICTIONARY(G__BCSUtil TChannel.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/Util/LinkDef.h)   
add_library(BCSUtil SHARED G__BCSUtil.cxx libraries/Util/TChannel.cxx)               
target_link_libraries(BCSUtil PUBLIC ${ROOT_LIBRARIES})                                                                  


ROOT_GENERATE_DICTIONARY(G__DDAS DDASEvent.h DDASHit.h ddaschannel.h  LINKDEF ${CMAKE_BINARY_DIR}/libraries/DDAS/LinkDef.h)   
FILE(GLOB DDASSources libraries/DDAS/*.cxx)
add_library(DDAS SHARED G__DDAS.cxx ${DDASSources})               
target_link_libraries(DDAS PUBLIC ${ROOT_LIBRARIES})                                                                  

ROOT_GENERATE_DICTIONARY(G__JUNK junk.h LINKDEF ${CMAKE_BINARY_DIR}/libraries/junk/LinkDef.h)   
add_library(Junk SHARED G__JUNK.cxx libraries/junk/junk.cxx)               
target_link_libraries(Junk PUBLIC ${ROOT_LIBRARIES})                                                                  


add_executable(bcssort src/bcssort.cxx) 

target_link_libraries(bcssort ${ROOT_LIBRARIES} BCSint BCSUtil DDAS Junk)



