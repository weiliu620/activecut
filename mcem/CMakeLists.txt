cmake_minimum_required(VERSION 2.8)
PROJECT(mcem)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

#add definitions, compiler switches, etc.
#ADD_DEFINITIONS(-O0)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_49_0)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/eigen-eigen-3.0.3)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})

ADD_EXECUTABLE(createimage 
  createimage.cxx
  util.cxx
  )

ADD_EXECUTABLE(createpotts
  createpotts.cxx
  util.cxx
  )

ADD_EXECUTABLE(segimage 
  segimage.cxx
  descent.cxx
  util.cxx
  )

ADD_EXECUTABLE(segimagealt 
  segimagealt.cxx
  MCModel.cxx
  GCoptimization.cpp 
  graph.cpp
  LinkedBlockList.cpp
  maxflow.cpp
  utilalt.cxx
  )

ADD_EXECUTABLE(verifyimage 
  verify.cxx
  )

ADD_EXECUTABLE(mlebeta
  mlebeta.cxx
  MCModel.cxx
  GCoptimization.cpp 
  graph.cpp
  LinkedBlockList.cpp
  maxflow.cpp
  utilalt.cxx
  )

ADD_EXECUTABLE(students_t
students_t_single_sample.cpp
  )


#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(createimage ITKIO ITKCommon )
TARGET_LINK_LIBRARIES(createpotts ITKIO ITKCommon )
TARGET_LINK_LIBRARIES(segimage ITKIO ITKCommon)
TARGET_LINK_LIBRARIES(segimagealt ITKIO ITKCommon )
TARGET_LINK_LIBRARIES(verifyimage ITKIO ITKCommon )
TARGET_LINK_LIBRARIES(mlebeta ITKIO ITKCommon)

