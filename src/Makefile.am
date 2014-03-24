#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS = validation

#------------------------------------
# Problem sources
#------------------------------------
validation_SOURCES = $(top_srcdir)/src/compute.C \
                     $(top_srcdir)/src/main.C 

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(QUESO_CPPFLAGS) \
               $(BOOST_CPPFLAGS) \
               $(GSL_CFLAGS) \
               $(HDF5_CFLAGS) \
               $(GLPK_CFLAGS) 

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(MPILIBS) \
             $(QUESO_LIBS) \
             $(GSL_LIBS) \
             $(GLPK_LIBS) \
             $(HDF5_LIBS) \
             $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
	           $(BOOST_PROGRAM_OPTIONS_LIBS) 
