configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape71"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape71" COPYONLY )
	     
configure_file("${RESOURCES}/n-040_Zr_090-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )
configure_file("${RESOURCES}/tsl-ZrinZrH-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape30" COPYONLY )

add_test( NAME "Test49"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
