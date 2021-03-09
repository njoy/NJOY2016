configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape31"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape31" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape32"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape32" COPYONLY )

configure_file("${RESOURCES}/n-026_Fe_000-IRDFF-II.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test60"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
