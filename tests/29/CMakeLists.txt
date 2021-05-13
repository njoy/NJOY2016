configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape35"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape35" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape45"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape45" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape55"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape55" COPYONLY )

configure_file("${RESOURCES}/n-094_Pu_241-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test29"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
