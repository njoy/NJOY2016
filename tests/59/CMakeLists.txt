configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape27"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape27" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape28"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape28" COPYONLY )

configure_file("${RESOURCES}/photoat-092_U_000-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

configure_file("${RESOURCES}/atrelax-092_U_000-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape21" COPYONLY )

configure_file("${RESOURCES}/photoat-092_U_000-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape40" COPYONLY )

configure_file("${RESOURCES}/atrelax-092_U_000-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape41" COPYONLY )

add_test( NAME "Test59"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
