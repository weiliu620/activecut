
PROJECT(gentoy)

SET(CMAKE_VERBOSE_MAKEFILE ON)

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -O0)

include_directories(/scratch/weiliu/rpms/boost_1_38_0)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

#list all source files here
ADD_EXECUTABLE(gentoy main.cpp GetNormalSample.cpp)

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(gentoy ITKIO ITKCommon)
