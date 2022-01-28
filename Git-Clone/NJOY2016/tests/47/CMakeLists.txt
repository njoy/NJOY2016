configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape45"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape45" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape46"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape46" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape47"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape47" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape48"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape48" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape55"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape55" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape56"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape56" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape57"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape57" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape58"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape58" COPYONLY )

configure_file("${RESOURCES}/n-094_Pu_239-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test47"
          COMMAND ${PYTHON_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
