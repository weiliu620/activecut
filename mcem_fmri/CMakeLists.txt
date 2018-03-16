cmake_minimum_required(VERSION 2.8)
PROJECT(mcem_fmri)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

#add definitions, compiler switches, etc.
#ADD_DEFINITIONS(-O0)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_46_1/build/include)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

LINK_DIRECTORIES (/home/sci/weiliu/packages/boost_1_46_1/build/lib)

MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})

ADD_EXECUTABLE(createmask 
  createmask.cxx
  utilalt.cxx
  )

ADD_EXECUTABLE(createfmri 
  createfmri.cxx
  vonmise.cxx
  utilalt.cxx
  )

ADD_EXECUTABLE(create3ddata
  create3ddata.cxx
  vonmise.cxx
  utilalt.cxx
  )

ADD_EXECUTABLE(projectfmri
  projectfmri.cxx
  utilalt.cxx
  )

ADD_EXECUTABLE(projecttoy
  projecttoy.cxx
  utilalt.cxx
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

ADD_EXECUTABLE(threshlabel 
  threshlabel.cxx
  )

ADD_EXECUTABLE(zshiftica
  zshiftica.cxx
  )

ADD_EXECUTABLE(getvar
  getvar.cxx
  )

ADD_EXECUTABLE(getsum
  getsum.cxx
  )

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(createmask ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(createfmri ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(create3ddata ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(segimagealt ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(projectfmri ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(projecttoy ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(threshlabel ${ITK_LIBRARIES} boost_program_options)
TARGET_LINK_LIBRARIES(zshiftica ${ITK_LIBRARIES} boost_program_options )
TARGET_LINK_LIBRARIES(getvar ${ITK_LIBRARIES} boost_program_options boost_filesystem)
TARGET_LINK_LIBRARIES(getsum ${ITK_LIBRARIES} boost_program_options)
