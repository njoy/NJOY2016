configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape50"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape50" COPYONLY )

configure_file("${RESOURCES}/n-018_Ar_37-tendl2023.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test85"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
