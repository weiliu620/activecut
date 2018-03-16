cmake_minimum_required(VERSION 2.8)
PROJECT(actlearn)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

option (USE_OMP
  "use openmp for parallel execution" on)
if (USE_OMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif(USE_OMP)


include_directories(${PROJECT_SOURCE_DIR})
set(Boost_DEBUG 0)
set(Boost_USE_STATIC_LIBS on)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options filesystem system)

if(Boost_FOUND)
  MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(grabcutseg
    grabcutseg.cxx
    loadfiles.cxx
    gmm.cxx
    logistic.cxx
    graphcuts.cxx
    graph.cpp
    maxflow.cpp
    utility.cxx
    )
  
  add_executable(copyheader
    copyheader.cxx
    )

  add_executable(pca
    pca.cxx
    utility.cxx
    )

  add_executable(rescale
    rescale.cxx
    )

  add_executable(queryccp
    queryccp.cxx
    query.cxx
    gmm.cxx
    loadfiles.cxx)

  target_link_libraries(grabcutseg ${ITK_LIBRARIES} ${Boost_LIBRARIES})
  target_link_libraries(copyheader ${ITK_LIBRARIES} ${Boost_LIBRARIES})
  target_link_libraries(pca ${ITK_LIBRARIES} ${Boost_LIBRARIES})
  target_link_libraries(rescale ${ITK_LIBRARIES} ${Boost_LIBRARIES})
  target_link_libraries(queryccp ${ITK_LIBRARIES} ${Boost_LIBRARIES})

endif()

