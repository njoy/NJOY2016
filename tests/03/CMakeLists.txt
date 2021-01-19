configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape37"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape37" COPYONLY )

configure_file("${RESOURCES}/gam23"
               "${CMAKE_CURRENT_BINARY_DIR}/tape30" COPYONLY )

configure_file("${RESOURCES}/gam27"
               "${CMAKE_CURRENT_BINARY_DIR}/tape32" COPYONLY )

add_test( NAME "Test03"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
