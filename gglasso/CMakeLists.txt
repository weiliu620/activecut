cmake_minimum_required(VERSION 2.8)
PROJECT(gglasso)

#add definitions, compiler switches, etc.

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)



# INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_43_0)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/eigen-eigen-3.0.3/)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_46_1/build/include)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

LINK_DIRECTORIES (/home/sci/weiliu/packages/boost_1_46_1/build/lib)


ADD_EXECUTABLE(testlasso
  testlasso.cxx
  codescent.cxx)

#need to link to some other libraries ? just add them here
TARGET_LINK_LIBRARIES(testlasso ${ITK_LIBRARIES} boost_program_options boost_filesystem)
