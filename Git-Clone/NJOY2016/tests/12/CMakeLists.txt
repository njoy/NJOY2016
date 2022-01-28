configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape22"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape22" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape24"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape24" COPYONLY )
	     
configure_file("${RESOURCES}/eni61"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test12"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
