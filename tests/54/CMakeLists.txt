configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape33"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape33" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape34"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape34" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape35"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape35" COPYONLY )
	     
configure_file("${RESOURCES}/p-001_H_003-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test54"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
