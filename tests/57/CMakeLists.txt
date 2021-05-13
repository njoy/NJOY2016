configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape31"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape31" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape34"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape34" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape37"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape37" COPYONLY )

configure_file("${RESOURCES}/g-083_Bi_209-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test57"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
