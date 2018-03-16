cmake_minimum_required(VERSION 2.8)
PROJECT(actlearn)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

include_directories(${PROJECT_SOURCE_DIR})

set(Boost_DEBUG 0)
set(Boost_USE_STATIC_LIBS on)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME OFF)
site_name(hostname)
set(Boost_INCLUDE_DIR /home/sci/weiliu/packages/boost_1_53_0/build/${hostname}/include)
set(Boost_LIBRARY_DIR /home/sci/weiliu/packages/boost_1_53_0/build/${hostname}/lib)
find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options filesystem system)

if(Boost_FOUND)
  MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(hdmrf
    hdmrf.cxx
    utility.cxx
    )
  target_link_libraries(hdmrf ${ITK_LIBRARIES} ${Boost_LIBRARIES})
endif()

