configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape27"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape27" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape28"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape28" COPYONLY )

configure_file("${RESOURCES}/t404"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test11"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
