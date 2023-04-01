configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape26"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape26" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape91"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape91" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape92"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape92" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape93"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape93" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape99"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape99" COPYONLY )

configure_file("${RESOURCES}/J33U238"
               "${CMAKE_CURRENT_BINARY_DIR}/tape21" COPYONLY )

configure_file("${RESOURCES}/J33U235"
               "${CMAKE_CURRENT_BINARY_DIR}/tape22" COPYONLY )

configure_file("${RESOURCES}/J33Pu239"
               "${CMAKE_CURRENT_BINARY_DIR}/tape23" COPYONLY )

add_test( NAME "Test17"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
