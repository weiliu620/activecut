cmake_minimum_required(VERSION 2.8)
PROJECT(utilities)

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

include_directories(${Boost_INCLUDE_DIR})

if(Boost_FOUND)
  ADD_LIBRARY(util
    util.cxx)

  ADD_EXECUTABLE(rain
    rain.cxx)	

  ADD_EXECUTABLE(compute_hist
    compute_hist.cxx)	

  ADD_EXECUTABLE(overlay
    overlay.cxx)	

  ADD_EXECUTABLE(overlayheat
    overlayheat.cxx)	

  ADD_EXECUTABLE(label_contour_filter
    label_contour_filter.cxx
    util.cxx)	

  ADD_EXECUTABLE(antialias_filter
    antialias_filter.cxx
    util.cxx)	

  TARGET_LINK_LIBRARIES(util  ${ITK_LIBRARIES} )
  TARGET_LINK_LIBRARIES(rain  ${Boost_LIBRARIES} ${ITK_LIBRARIES} )
  TARGET_LINK_LIBRARIES(compute_hist  ${Boost_LIBRARIES}  ${ITK_LIBRARIES} )
  TARGET_LINK_LIBRARIES(overlay  ${Boost_LIBRARIES}  ${ITK_LIBRARIES} )
  TARGET_LINK_LIBRARIES(overlayheat  ${Boost_LIBRARIES}  ${ITK_LIBRARIES} )
  TARGET_LINK_LIBRARIES(label_contour_filter ${Boost_LIBRARIES}  ${ITK_LIBRARIES} )
endif()






