#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(PredictiveSelectionCoupledModels, 0.1.0, terejanu@ices.utexas.edu)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE(PredictiveSelectionCoupledModels,0.1.0)

AC_PROG_LIBTOOL

#------------------------------
# Checks for compilers and MPI
#------------------------------
AC_PROG_CC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_PROG_CXX
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

#---------------------------------------------------------
# Checks for library prerequisites for other libraries...
#---------------------------------------------------------

#----------------
# Check for GSL
#----------------
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))

#----------------
# Check for boost
#----------------
BOOST_REQUIRE([1.36])
BOOST_MATH()
BOOST_PROGRAM_OPTIONS()

#----------------
# Check for GLPK (optional)
#----------------
AX_PATH_GLPK([4.35],[no])

#----------------
# Check for HDF5 (optional)
#----------------
AX_PATH_HDF5_NEW([1.8.0],[no])

#----------------
# Check for QUESO with ANN
#----------------
AX_PATH_QUESO_ANN

#-------------------
# Generate makefiles
#-------------------
AC_OUTPUT( 	Makefile \
	   		src/Makefile \
	   		tests/Makefile )

#--------------------------------
# Query configuration environment
#------------------------------
AX_SUMMARIZE_ENV

#------------------------------
# Final summary
#------------------------------
AX_SUMMARIZE_CONFIG

