# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: src/test/CMakeLists.txt                               #
#                                                             #
# Project name: GiNaCRA - GiNaC Real Algebra package          #
# Contact: Ulrich Loup <loup@cs.rwth-aachen.de>               #
# Project sourceforge page: http://ginacra.sourceforge.net/   #
# Development status: beta (testing + some features missing)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 2012-03-15                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include_directories( "${PROJECT_SOURCE_DIR}/src/lib" )
link_directories( "${PROJECT_BINARY_DIR}/src/lib" )

set( check_ginacra_SOURCES
     utilities_unittest.h utilities_unittest.cpp
     OpenInterval_unittest.h OpenInterval_unittest.cpp
	 UnivariatePolynomial_unittest.h UnivariatePolynomial_unittest.cpp RationalUnivariatePolynomial_unittest.h RationalUnivariatePolynomial_unittest.cpp
	 RealAlgebraicNumberIR_unittest.h RealAlgebraicNumberIR_unittest.cpp RealAlgebraicNumber_unittest.h RealAlgebraicNumber_unittest.cpp RealAlgebraicNumberFactory_unittest.h RealAlgebraicNumberFactory_unittest.cpp
	 Constraint_unittest.h Constraint_unittest.cpp
	 CAD_unittest.h CAD_unittest.cpp
	 MultivariateMonomialMR_unittest.h MultivariateMonomialMR_unittest.cpp MultivariateTermMR_unittest.h MultivariateTermMR_unittest.cpp MultivariatePolynomialMR_unittest.h MultivariatePolynomialMR_unittest.cpp Groebner_unittest.h Groebner_unittest.cpp
	 test.cpp )

# main target
add_executable( check_${PROJECT_NAME}
                ${check_ginacra_SOURCES} )

# executable target
add_custom_target( check
                   COMMAND check_${PROJECT_NAME}
                   DEPENDS check_${PROJECT_NAME}
                   COMMENT "${ginacra_NAME} ${ginacra_VERSION} CPPUnit Test Suite" )

# link the libraries
target_link_libraries( check_${PROJECT_NAME}
                       ${CPPUNIT_LIBRARIES}
                       ${cln_LIBRARIES}
                       ${GINAC_LIBRARY}
                       ${PROJECT_NAME} )
