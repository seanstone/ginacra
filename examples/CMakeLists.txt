# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: examples/CMakeLists.txt                               #
#                                                             #
# Project name: GiNaCRA - GiNaC Real Algebra package          #
# Contact: Ulrich Loup <loup@cs.rwth-aachen.de>               #
# Project sourceforge page: http://ginacra.sourceforge.net/   #
# Development status: beta (testing + some features missing)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 2012-05-19                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include_directories( "${PROJECT_SOURCE_DIR}/src/lib" )
include_directories( "${PROJECT_BINARY_DIR}/include" ) # in order to find ginacra headers
include_directories( "${PROJECT_BINARY_DIR}" ) # in order to find config.h
link_directories( "${PROJECT_BINARY_DIR}/src/lib" )

# list all examples to be build
# set( EXAMPLE_DIR = ${PROJECT_SOURCE_DIR}/examples )
# set( EXAMPLES EXAMPLE_DIR/*.cpp )
set( EXAMPLES benchmark_Groebner_1
              example1 example2 example3 example4
              example_AlgoSyn-presentation-2012-07
#              example_NFM2011-presentation
              example_CAD_1 example_CAD_2 example_CAD_3 example_CAD_4 example_CAD_5 example_CAD_6 example_CAD_7 example_CAD_8 example_CAD_9 example_CAD_A example_CAD_B example_CAD_C
              example_CADincremental_1 example_CADincremental_2 example_CADincremental_3
              example_IR_1 example_IR_2 example_IR_3 example_IR_4
              example_RAN_1 example_RAN_2 example_RAN_3
               )

set( EXAMPLE_SOURCES benchmark_Groebner_1.cpp
                     example1.cpp example2.cpp example3.cpp example4.cpp
                     example_AlgoSyn-presentation-2012-07.cpp
#                     example_NFM2011-presentation.cpp
                     example_CAD_1.cpp example_CAD_2.cpp example_CAD_3.cpp example_CAD_4.cpp example_CAD_5.cpp example_CAD_6.cpp example_CAD_7.cpp example_CAD_8.cpp example_CAD_9.cpp  example_CAD_A.cpp example_CAD_B.cpp example_CAD_C.cpp
                     example_CADincremental_1 example_CADincremental_2 example_CADincremental_3
                     example_IR_1.cpp example_IR_2.cpp example_IR_3.cpp example_IR_4.cpp
                     example_RAN_1.cpp example_RAN_2.cpp
#                     example_SGB_1
               )

foreach( example ${EXAMPLES} )
    add_executable( ${example} ${example}.cpp )
    target_link_libraries( ${example} ${cln_LIBRARIES} ${GINAC_LIBRARY} ${PROJECT_NAME} ${READLINE_LIBRARY} )
endforeach( example )

add_custom_target( examples
                   DEPENDS headers
                   DEPENDS ${EXAMPLES}
                   COMMENT "${ginacra_NAME} ${ginacra_VERSION} Examples" )
