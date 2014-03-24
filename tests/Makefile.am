#------------------------------------
# Tests we want to build
#------------------------------------
check_PROGRAMS = 	forcing_SLD_funcs_test          \
		 			forcing_SED_funcs_test          \
		 			forcing_OLD_funcs_test          \
		 			forcing_OED_funcs_test          \
	         		oscillator_OLS_funcs_test       \
	         		oscillator_OCS_funcs_test       \
                 	oscillator_OQS_funcs_test       \
                 	coupled_OLS_SED_funcs_test

#------------------------------------
# Sources for these tests
#------------------------------------
forcing_SLD_funcs_test_SOURCES = forcing_SLD_funcs_test.C
forcing_SED_funcs_test_SOURCES = forcing_SED_funcs_test.C
forcing_OLD_funcs_test_SOURCES = forcing_OLD_funcs_test.C
forcing_OED_funcs_test_SOURCES = forcing_OED_funcs_test.C
oscillator_OLS_funcs_test_SOURCES = oscillator_OLS_funcs_test.C
oscillator_OCS_funcs_test_SOURCES = oscillator_OCS_funcs_test.C
oscillator_OQS_funcs_test_SOURCES = oscillator_OQS_funcs_test.C
coupled_OLS_SED_funcs_test = coupled_OLS_SED_funcs_test.C

#------------------------------------
# Define tests to actually be run
#------------------------------------
TESTS = 	models.sh                    \
			forcing_SLD_funcs_test       \
			forcing_SED_funcs_test       \
			forcing_OLD_funcs_test       \
			forcing_OED_funcs_test       \
			oscillator_OLS_funcs_test    \
			oscillator_OCS_funcs_test    \
	        oscillator_OQS_funcs_test    \
	        coupled_OLS_SED_funcs_test

#------------------------------------
# Want these put with the distro so we can run make check
#------------------------------------
shellfiles_src = models.sh

EXTRA_DIST = $(shellfiles_src)

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include   \
	           -I$(top_srcdir)/tests

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(HDF5_CFLAGS) \
               $(QUESO_CPPFLAGS) \
               $(BOOST_CPPFLAGS) \
               $(GSL_CFLAGS) \
               $(GLPK_CFLAGS) 

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(MPILIBS) \
             $(QUESO_LIBS) \
             $(GSL_LIBS) \
             $(GLPK_LIBS) \
             $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
	         $(BOOST_PROGRAM_OPTIONS_LIBS) \
             $(HDF5_LIBS) 

