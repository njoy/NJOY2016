configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape60"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape60" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape55"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape55" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape71"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape71" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape76"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape76" COPYONLY )

configure_file("${RESOURCES}/n-001_H_002-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

configure_file("${RESOURCES}/tsl-DinLiD-ENDF8.1-beta.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape21" COPYONLY )

add_test( NAME "Test67"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
