# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: src/lib/CMakeLists.txt                                #
#                                                             #
# Project name: GiNaCRA - GiNaC Real Algebra package          #
# Contact: Ulrich Loup <loup@cs.rwth-aachen.de>               #
# Project sourceforge page: http://ginacra.sourceforge.net/   #
# Development status: beta (testing + some features missing)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 2012-03-20                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set( lib_core_headers
     ginacra.h tree.h constants.h utilities.h settings.h
     OpenInterval.h
     Polynomial.h UnivariatePolynomial.h UnivariatePolynomialSet.h RationalUnivariatePolynomial.h
     Constraint.h
     RealAlgebraicNumber.h RealAlgebraicNumberNR.h RealAlgebraicNumberIR.h operators.h RealAlgebraicNumberFactory.h RealAlgebraicPoint.h
     CAD.h )

set( lib_core_src
     ${lib_core_headers}
     utilities.cpp settings.cpp
     OpenInterval.cpp
     Polynomial.cpp UnivariatePolynomial.cpp UnivariatePolynomialSet.cpp RationalUnivariatePolynomial.cpp
     Constraint.cpp
     RealAlgebraicNumber.cpp RealAlgebraicNumberIR.cpp RealAlgebraicNumberNR.cpp operators.cpp RealAlgebraicNumberFactory.cpp RealAlgebraicPoint.cpp
     CAD.cpp )

set( lib_MRBranch_src
     SymbolDB.h VariableListPool.h
     SymbolDB.cpp VariableListPool.cpp
     MultivariateMonomialMR.h MultivariateTermMR.h MultivariatePolynomialMR.h MultivariateCoefficientMR.h Groebner.h
     MultivariateMonomialMR.cpp MultivariateTermMR.cpp MultivariatePolynomialMR.cpp MultivariateCoefficientMR.cpp Groebner.cpp )

set( lib_${PROJECT_NAME}_SRCS
     ${lib_core_src} ${lib_MRBranch_src} )

# determine whether we build a static library
set(STATICLIB "SHARED")
if(${STATICLIB_SWITCH} STREQUAL "ON")
	set(STATICLIB "STATIC")
endif()
add_library( lib_${PROJECT_NAME}
             ${STATICLIB}
             ${lib_${PROJECT_NAME}_SRCS} )

# target to install the headers into the build directory
add_custom_target( headers install -d ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME} && install ${lib_core_headers} ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}/
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/lib
                   COMMENT "Copy headers to build directory" )

# other libraries that are used by the project library should go here.
target_link_libraries( lib_${PROJECT_NAME} ${cln_LIBRARIES} ${GINAC_LIBRARY} )

set_target_properties( lib_${PROJECT_NAME} PROPERTIES
#		# create *nix style library versions + symbolic links
		VERSION "${ginacra_MAJORVERSION}.${ginacra_MINORVERSION}.${ginacra_MAINTENANCEVERSION}"
		SOVERSION "${ginacra_MAJORVERSION}.${ginacra_MINORVERSION}.${ginacra_MAINTENANCEVERSION}"
#		# allow creating static and shared libs without conflicts
#		# CLEAN_DIRECT_OUTPUT 1
#		# avoid conflicts between library and binary target names
		OUTPUT_NAME ${PROJECT_NAME} )

# set install destination and permissions for the ginacra library
install( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
         OWNER_READ OWNER_WRITE OWNER_EXECUTE
         GROUP_READ GROUP_EXECUTE
         WORLD_READ WORLD_EXECUTE )

install( FILES
         ginacra.h tree.h constants.h utilities.h settings.h
		 OpenInterval.h
		 Polynomial.h UnivariatePolynomial.h UnivariatePolynomialSet.h RationalUnivariatePolynomial.h
		 Constraint.h
		 RealAlgebraicNumber.h RealAlgebraicNumberNR.h RealAlgebraicNumberIR.h operators.h RealAlgebraicNumberFactory.h
		 RealAlgebraicPoint.h
		 CAD.h
		 SymbolDB.h VariableListPool.h
		 MultivariateMonomialMR.h MultivariateTermMR.h MultivariatePolynomialMR.h MultivariateCoefficientMR.h Groebner.h
         DESTINATION include/ginacra )
