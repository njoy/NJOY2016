configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape31"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape31" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape34"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape34" COPYONLY )

configure_file("${RESOURCES}/g-002-He-3-IAEA.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test78"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
