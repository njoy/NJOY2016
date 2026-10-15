configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape30"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape30" COPYONLY )

configure_file("${RESOURCES}/n-038_Sr_088-ENDF8.1.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test81"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
