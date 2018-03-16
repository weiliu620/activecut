cmake_minimum_required(VERSION 2.8)
PROJECT(graphcuts)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-O0)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

set(Boost_DEBUG 0)
set(Boost_USE_STATIC_LIBS on)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME OFF)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/eigen-eigen-ffa86ffb5570)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/gco/)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

site_name(hostname)
set(Boost_INCLUDE_DIR /home/sci/weiliu/packages/boost_1_53_0/build/${hostname}/include)
set(Boost_LIBRARY_DIR /home/sci/weiliu/packages/boost_1_53_0/build/${hostname}/lib)
find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options filesystem system)

include_directories(${Boost_INCLUDE_DIRS})


LINK_DIRECTORIES (/home/sci/weiliu/packages/boost_1_43_0/build/lib)

ADD_EXECUTABLE(generateimage 
  generateimage.cxx
  util.cxx)

ADD_EXECUTABLE(recoverimage 
  recoverimage.cxx 
  GCoptimization.cpp 
  graph.cpp
  LinkedBlockList.cpp
  maxflow.cpp
  util.cxx)

# ADD_LIBRARY(graphcuts
#   GCoptimization.cpp 
#   graph.cpp
#   LinkedBlockList.cpp
#   maxflow.cpp
#   util.cxx)

ADD_EXECUTABLE(compenergy 
  compenergy.cxx
  GCoptimization.cpp 
  graph.cpp
  LinkedBlockList.cpp
  maxflow.cpp
  util.cxx)

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(generateimage ${ITK_LIBRARIES} ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(recoverimage ${ITK_LIBRARIES} ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(compenergy ${ITK_LIBRARIES} ${Boost_LIBRARIES} ) 
# TARGET_LINK_LIBRARIES(graphcuts ${ITK_LIBRARIES}  ${Boost_LIBRARIES} )
