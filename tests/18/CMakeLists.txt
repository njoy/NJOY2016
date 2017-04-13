configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape28"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape28" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape38"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape38" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape39"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape39" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape91"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape91" COPYONLY )
	     
configure_file("${RESOURCES}/DCf252"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test18"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
