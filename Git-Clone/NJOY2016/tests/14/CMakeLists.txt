configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape31"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape31" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape36"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape36" COPYONLY )

configure_file("${RESOURCES}/epn14"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

configure_file("${RESOURCES}/epn14"
               "${CMAKE_CURRENT_BINARY_DIR}/tape21" COPYONLY )

add_test( NAME "Test14"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
