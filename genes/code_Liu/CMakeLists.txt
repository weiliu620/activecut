cmake_minimum_required(VERSION 2.8)
PROJECT(groupmrf)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_53_0)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.2.3)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon_build)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

LINK_DIRECTORIES (/home/sci/weiliu/packages/lemon_build/lemon)

# boost library
set(BOOST_LIB "" CACHE FILEPATH "Boost library path")
LINK_DIRECTORIES( ${BOOST_LIB} )

ADD_EXECUTABLE(geteigenimage
  geteigenimage.cxx
  )


TARGET_LINK_LIBRARIES(geteigenimage ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
