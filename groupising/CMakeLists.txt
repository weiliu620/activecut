cmake_minimum_required(VERSION 2.8)
PROJECT(groupising)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

#add definitions, compiler switches, etc.

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_53_0)

INCLUDE_DIRECTORIES (/home/sci/weiliu/projects/util)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
include_directories ("${PROJECT_SOURCE_DIR}/gco")

# boost library
set(BOOST_LIB "" CACHE FILEPATH "Boost library path")
LINK_DIRECTORIES( ${BOOST_LIB} )


LINK_DIRECTORIES (/home/sci/weiliu/projects/util)
LINK_DIRECTORIES ("${PROJECT_SOURCE_DIR}/gco")

MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})


ADD_EXECUTABLE(createmask
  createmask.cxx
  util.cxx
)

ADD_EXECUTABLE(creategrpmap
  createimage.cxx
  util.cxx
  )
ADD_EXECUTABLE(gensubmap
  gensubmap.cxx
  util.cxx
  )

ADD_EXECUTABLE(validatelabelmap
  validatelabelmap.cxx
  util.cxx
  )

ADD_EXECUTABLE(genobsimage
  genobsimage.cxx
  )

ADD_EXECUTABLE(gmmseg
  gmmseg.cxx
  gmmodel.cxx
  util.cxx
  )

ADD_EXECUTABLE(kmeans
  kmeans.cxx
  util.cxx
  )

# ADD_EXECUTABLE(mcemgmm
#   mcemgmm.cxx
#   gmmodel_mcem.cxx
#   gmmodel_graph.cxx
#   util.cxx
#   )

# ADD_EXECUTABLE(parest
#   parest.cxx
#   )

ADD_EXECUTABLE(genest
  genest.cxx
  util.cxx
  )

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(createmask ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
TARGET_LINK_LIBRARIES(creategrpmap ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
TARGET_LINK_LIBRARIES(gensubmap ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
TARGET_LINK_LIBRARIES(validatelabelmap ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
TARGET_LINK_LIBRARIES(genobsimage ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
TARGET_LINK_LIBRARIES(gmmseg ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
TARGET_LINK_LIBRARIES(kmeans ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)
# TARGET_LINK_LIBRARIES(mcemgmm ${ITK_LIBRARIES} graphcuts boost_program_options boost_filesystem boost_system)

# TARGET_LINK_LIBRARIES(parest ${ITK_LIBRARIES} graphcuts boost_program_options boost_filesystem boost_system)

TARGET_LINK_LIBRARIES(genest ${ITK_LIBRARIES}  boost_program_options boost_filesystem boost_system)


# add_subdirectory (gco) 

SUBDIRS(${groupising_SOURCE_DIR}/undigraph)



