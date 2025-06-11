configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape25"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape25" COPYONLY )

configure_file("${RESOURCES}/n-Hf177-TENDL2021.tendl"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )

add_test( NAME "Test86"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
