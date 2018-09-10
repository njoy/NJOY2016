configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape25"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape25" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape26"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape26" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape27"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape27" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape35"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape35" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape36"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape36" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape37"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape37" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape45"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape45" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape46"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape46" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape47"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape47" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape91"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape91" COPYONLY )
	     
configure_file("${RESOURCES}/J33U238"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test15"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
