configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape43"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape43" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape44"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape44" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape45"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape45" COPYONLY )

configure_file("${RESOURCES}/n-001_H_001-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )
configure_file("${RESOURCES}/photoat-001_H_000-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape40" COPYONLY )

add_test( NAME "Test30"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
