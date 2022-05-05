configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape41"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape41" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape44"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape44" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape54"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape54" COPYONLY )

configure_file("${RESOURCES}/n-004_Be_009-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test72"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
