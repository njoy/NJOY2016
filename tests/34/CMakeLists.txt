configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape33"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape33" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape43"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape43" COPYONLY )
	     
configure_file("${RESOURCES}/n-094_Pu_240-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test34"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
