configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape71"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape71" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape72"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape72" COPYONLY )

configure_file("${RESOURCES}/thermalscattering.ace"
               "${CMAKE_CURRENT_BINARY_DIR}/tape61" COPYONLY )

add_test( NAME "Test61"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
