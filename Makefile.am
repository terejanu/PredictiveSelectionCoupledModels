AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS  = -I m4

SUBDIRS          = src tests

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES = aclocal.m4                               \
                       aminclude.am                             \
                       autom4te.cache/*                         \
                       $(AUX_DIST)                              \
                       config.log                               \
                       config.status                            \
                       config.sub                               \
                       configure                                \
                       config.h                                 \
                       config.h.in                              \
                       Makefile.in                              \
                       src/Makefile.in                          \
		                   tests/Makefile.in

# Variables for grouping the experimental data and
# input files. 

input_file_SRC = $(top_srcdir)/input_files/queso.inp \
	               $(top_srcdir)/input_files/forcing_meas.txt \
                 $(top_srcdir)/input_files/oscillator_meas.txt

matlab_testConv_SRC = $(top_srcdir)/matlab/*.m
matlab_genData_SRC = $(top_srcdir)/matlab/generateData/*.m
matlab_mcmcdiag_SRC = $(top_srcdir)/matlab/mcmcdiag/*.*

# Want the experimental data and input files to be distributed
# We add the regession test experimental data for running
# make check on a distribution, but we don't need them
# installed.
EXTRA_DIST = $(input_file_SRC)

# Want the experimental data and input files to be installed

input_filedir = $(prefix)/input_files
input_file_DATA = $(input_file_SRC)

matlab_testConvdir = $(prefix)/matlab
matlab_testConv_DATA = $(matlab_testConv_SRC)

matlab_genDatadir = $(prefix)/matlab/generateData
matlab_genData_DATA = $(matlab_genData_SRC)

matlab_mcmcdiagdir = $(prefix)/matlab/mcmcdiag
matlab_mcmcdiag_DATA = $(matlab_mcmcdiag_SRC)
