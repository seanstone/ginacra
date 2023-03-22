# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: src/tools/CMakeLists.txt                              #
#                                                             #
# Project name: GiNaCRA - GiNaC Real Algebra package          #
# Contact: Ulrich Loup <loup@cs.rwth-aachen.de>               #
# Project sourceforge page: http://ginacra.sourceforge.net/   #
# Development status: beta (testing + some features missing)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 2012-02-22                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include_directories( "${PROJECT_SOURCE_DIR}/src/lib" )
include_directories("${PROJECT_BINARY_DIR}") # in order to find config.h
link_directories( "${PROJECT_BINARY_DIR}/src/lib" )

set( ${PROJECT_NAME}console_SOURCES ${PROJECT_NAME}console.cpp )

add_executable( ${PROJECT_NAME}console EXCLUDE_FROM_ALL ${${PROJECT_NAME}console_SOURCES} )

# Link the libraries
target_link_libraries( ${PROJECT_NAME}console
                       ${cln_LIBRARIES}
                       ${GINAC_LIBRARY}
                       ${PROJECT_NAME}
                       ${READLINE_LIBRARY} )
