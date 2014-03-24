Predictive Selection of Coupled Models
=======================

This is an application framework for solving the problem of
predictive model selection of coupled models as presented in the 
following paper:

Gabriel Terejanu, Todd Oliver, Chris Simmons (2011). "Application of 
Predictive Model Selection to Coupled Models". In Proceedings of the World 
Congress on Engineering and Computer Science 2011 Vol II, WCECS 2011, 
pp. 927-932.

The framework is built on top of the statistical library QUESO 
(Quantification of Uncertainty for Estimation, Simulation and Optimization):
https://red.ices.utexas.edu/projects/software/wiki/QUESO

Dependencies
=======================

This software is dependent on QUESO 0.46 or later, which needs to 
be compiled along with ANN - a library for Approximate Nearest Neighbor 
Searching:
http://www.cs.umd.edu/~mount/ANN/

QUESO 0.46 contains a copy of ANN and can be configured using the
option --enable-ann=yes.

Other dependencies are similar with QUESO's dependencies, namely:
MPI, GSL, BOOST, GLPK, and HDF5.

Installation
=======================

This package uses GNU Autoconf system for configuration. Please provide
the installation location $INSTALL_DIR, and the location where QUESO
is installed, $QUESO_DIR.

$ ./bootstrap

$ ./configure --prefix=$INSTALL_DIR --with-queso=$QUESO_DIR

To build, check and install the software type the following commands.
Note that "make check" is optional. 

$ make

$ make check

$ make install 