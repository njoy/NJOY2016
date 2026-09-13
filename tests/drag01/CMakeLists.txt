configure_file("${CMAKE_CURRENT_SOURCE_DIR}/input"
               "${CMAKE_CURRENT_BINARY_DIR}/input" COPYONLY )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/referenceTape30"
               "${CMAKE_CURRENT_BINARY_DIR}/referenceTape30" COPYONLY )

configure_file("${RESOURCES}/tsl-HinH2O-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape40" COPYONLY )

configure_file("${RESOURCES}/n-001_H_001-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape41" COPYONLY )

configure_file("${RESOURCES}/n-026_Fe_058-ENDF8.0-Beta4.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape43" COPYONLY )

configure_file("${RESOURCES}/n-094_Pu_239-ENDF8.0-Beta6.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape44" COPYONLY )

configure_file("${RESOURCES}/n-036_Kr_083-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape46" COPYONLY )

configure_file("${RESOURCES}/n-050_Sn_119-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape47" COPYONLY )

configure_file("${RESOURCES}/nfy-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape50" COPYONLY )

configure_file("${RESOURCES}/decay-ENDF8.0.endf"
               "${CMAKE_CURRENT_BINARY_DIR}/tape51" COPYONLY )

add_test( NAME "drag01"
          COMMAND ${Python3_EXECUTABLE} "../execute.py"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
