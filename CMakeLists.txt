# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: CMakeLists.txt                                        #
#                                                             #
# Project name: GiNaCRA - GiNaC Real Algebra package          #
# Contact: Ulrich Loup <loup@cs.rwth-aachen.de>               #
# Project sourceforge page: http://ginacra.sourceforge.net/   #
# Development status: beta (testing + some features missing)  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 2012-05-19                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# common settings
# # # # # # # # #

cmake_minimum_required( VERSION 2.8 )
project( ginacra )

#set( CMAKE_BUILD_TYPE Release)
#set( CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -std=gnu++0x")

# set the install path
if( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "install prefix" FORCE )
endif( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )

# information about the project
set( PROJECT_FULLNAME "GiNaCRA")
set( PROJECT_DESCRIPTION "GiNaC Real Algebra Package")
set( ginacra_MAJORVERSION 0)
set( ginacra_MINORVERSION 6)
set( ginacra_MAINTENANCEVERSION 4)
set( ginacra_VERSION "${ginacra_MAJORVERSION}.${ginacra_MINORVERSION}.${ginacra_MAINTENANCEVERSION}")
set( PROJECT_VERSION ${ginacra_VERSION} )
set( ginacra_NAME \"${PROJECT_NAME}\")

# options
option( STATICLIB_SWITCH "Compile a statically linked version of the library." )
# if(${STATICLIB_SWITCH} STREQUAL "ON")
# 	set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static" )
# endif()


# path to find own modules
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )


# required libraries
# # # # # # # # # # #

# GiNaCRA libraries and tools all need GiNaC
# TODO: add support for version number
find_package( GiNaC 1.5.8 REQUIRED )

# find CLN ..
# TODO: we should write a find-module for this
include( FindPkgConfig )
pkg_check_modules( cln REQUIRED cln>=1.2.2 )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${cln_LDFLAGS}" )
# On ELF systems (such as GNU/Linux) this prevents the infamous `error
# while loading shared libraries: libfoo.so: cannot open shared object
# file: no such file or directory'
set( CMAKE_INSTALL_RPATH_USE_LINK_RPATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${cln_LIBRARY_DIRS}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
	set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${cln_LIBRARY_DIRS}")
endif()


# optional libraries
# # # # # # # # # # #

# find CPPUnit (optional), providing CppUnit_FOUND
find_package( CPPUnit )
# find readline (optional), providing Readline_FOUND
include( FindReadline )


# further directories to be included
# # # # # # # # # # # # # # # # # # #

add_subdirectory(src)
add_subdirectory(examples EXCLUDE_FROM_ALL) # offer as an optional target

# If doxygen is installed, create the make doc target.
include(FindDoxygen)
if(DOXYGEN)
  add_subdirectory(doxygen-conf)
else()
  message(WARNING, "Doxygen not found")
endif()


# configure file specification in order to translate CMake variables to C++
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set( DOLLAR "$" )
configure_file( ${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_BINARY_DIR}/config.h )
configure_file( ${CMAKE_SOURCE_DIR}/ginacra.pc.in ${CMAKE_BINARY_DIR}/ginacra.pc )
include_directories( ${CMAKE_BINARY_DIR} )


# [un]install targets
# # # # # # # # # # #

# specify pkg-config file to be installed
install( FILES ${CMAKE_BINARY_DIR}/ginacra.pc DESTINATION lib/pkgconfig )

# uninstall target (see http://www.cmake.org/Wiki/RecipeAddUninstallTarget)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeUninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/CMakeUninstall.cmake"
  IMMEDIATE @ONLY )
add_custom_target(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/CMakeUninstall.cmake" )


# packaging target (cpack)
# # # # # # # # # # # # # #

# set( CPACK_SET_DESTDIR ON )
set( CPACK_GENERATOR "TGZ" )
# set( CPACK_PACKAGE_FILE_NAME ginacra-${ginacra_VERSION})
set( CPACK_PACKAGE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README" )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION} )
set( CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING" )
set( CPACK_PACKAGE_VERSION_MAJOR ${ginacra_MAJORVERSION})
set( CPACK_PACKAGE_VERSION_MINOR ${ginacra_MINORVERSION})
set( CPACK_PACKAGE_VERSION_PATCH ${ginacra_MAINTENANCEVERSION})
# set( CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)

set( CPACK_SOURCE_GENERATOR "TGZ;TBZ2" )
set( CPACK_SOURCE_PACKAGE_FILE_NAME ginacra-${ginacra_VERSION} )
# set( CPACK_SOURCE_STRIP_FILES "")
set( CPACK_SOURCE_IGNORE_FILES "/build/;/htdocs/;/nbproject/;/CVS/;/\\\\.svn/;/\\\\.bzr/;/\\\\.hg/;/\\\\.git/;/\\\\.deps/;/\\\\.libs/;\\\\.swp$;\\\\.pc$;\\\\.xjs$;\\\\.sh$;ylwrap;\\\\.prej$;\\\\.bak$;\\\\.backup$;~$" )
set( CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES};
     MultivariatePolynomial.h;MultivariatePolynomial.tpp;MultivariateTerm.h;MultivariateTerm.tpp;MultivariateMonomial.h;
     MultivariateMonomial.tpp;MultivariateCoefficient.h;MultivariateCoefficient.tpp;
	 MultivariatePolynomialFactory.h;SpecialQuotientRingMultiplicationTable.h;SpecialQuotientRingMultiplicationTable.tpp
   )
set( CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES};
     MultivariatePolynomial_unittest.h;MultivariatePolynomial_unittest.cpp;MultivariateTerm_unittest.h;MultivariateTerm_unittest.cpp;MultivariateMonomial_unittest.h;
     MultivariateMonomial_unittest.cpp;MultivariateCoefficient_unittest.h;MultivariateCoefficient_unittest.cpp;
	 MultivariatePolynomialFactory_unittest.h;MultivariatePolynomialFactory_unittest.cpp;SpecialQuotientRingMultiplicationTable_unittest.h;SpecialQuotientRingMultiplicationTable_unittest.tpp
   )
# message("CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)
