cmake_minimum_required(VERSION 2.8)
PROJECT(mrf)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-O0)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)
# Find Cuda. Dont have to include anything like itk?
FIND_PACKAGE(CUDA REQUIRED)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/ImageViewer)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_42_0)

MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})
MESSAGE(STATUS "CUDA_VERSION: " ${CUDA_VERSION})

# add library path.
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/ImageViewer)

# add library path.
#LINK_DIRECTORIES("/scratch/apps/itk_app_shared/Auxiliary/FltkImageViewer")
LINK_DIRECTORIES("/home/sci/weiliu/packages/itk_app_shared/Auxiliary/FltkImageViewer")

#list all source files here
CUDA_ADD_EXECUTABLE(base base.cu)
CUDA_ADD_EXECUTABLE(confield main.cxx 
  ConfigFile.cpp
  ExtractRealData.cxx 
  ExtractToyData.cxx 
  EMPost.cxx 
  GPUCorr.cu 
  GPUInit.cu 
  GPUSampling.cu 
  GPUEst.cu 
  GPUDetrend.cu
  GPUAnnealing.cu
  tools.cxx
  graymatter.cxx)

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(confield ITKIO ITKCommon FieldViewer)
ADD_SUBDIRECTORY(ImageViewer) 

