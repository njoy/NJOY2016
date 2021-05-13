configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape24"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape24" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape34"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape34" COPYONLY )

add_test( NAME "Test33"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
