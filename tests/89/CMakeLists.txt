configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape30"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape30" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape60"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape60" COPYONLY )

configure_file("${RESOURCES}/n-001_H_001-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )
configure_file("${RESOURCES}/n-018_Ar_37-tendl2023.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape50" COPYONLY )

add_test( NAME "Test89"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
