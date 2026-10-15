configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape26"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape26" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape27"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape27" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape28"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape28" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape29"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape29" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape36"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape36" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape37"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape37" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape38"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape38" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape39"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape39" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape46"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape46" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape47"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape47" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape48"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape48" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape49"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape49" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape56"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape56" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape57"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape57" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape58"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape58" COPYONLY )
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape59"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape59" COPYONLY )

configure_file("${RESOURCES}/n-027_Co_058-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape20" COPYONLY )
configure_file("${RESOURCES}/n-027_Co_058m1-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape30" COPYONLY )
configure_file("${RESOURCES}/n-095_Am_242-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape40" COPYONLY )
configure_file("${RESOURCES}/n-095_Am_242m1-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape50" COPYONLY )

add_test( NAME "Test82"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
