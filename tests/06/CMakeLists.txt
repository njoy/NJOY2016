configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape32"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape32" COPYONLY )

configure_file("${RESOURCES}/t511"
               "${CMAKE_CURRENT_BINARY_DIR}/tape30" COPYONLY )

add_test( NAME "Test06"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
