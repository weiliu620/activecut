ADD_EXECUTABLE(gentwolayers
  gentwolayers.cxx
  ../util.cxx
)

TARGET_LINK_LIBRARIES(gentwolayers util ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)