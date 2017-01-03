configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape27"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape27" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape29"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape29" COPYONLY )
	     
configure_file("${RESOURCES}/e6pu241c"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test19"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
