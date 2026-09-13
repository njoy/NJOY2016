configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape51"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape51" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape52"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape52" COPYONLY )

configure_file("${RESOURCES}/n-094_Pu_239-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape40" COPYONLY )

configure_file("${RESOURCES}/n-094_Pu_240-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape41" COPYONLY )

add_test( NAME "calendf01"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
