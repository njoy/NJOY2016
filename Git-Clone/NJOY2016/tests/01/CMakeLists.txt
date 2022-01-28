configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape25"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape25" COPYONLY )

configure_file("${RESOURCES}/t511"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

configure_file("${RESOURCES}/t322"
               "${CMAKE_CURRENT_BINARY_DIR}/tape26" COPYONLY )

add_test( NAME "Test01"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
