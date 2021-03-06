/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 Gabriel Terejanu - terejanu@cec.sc.edu
 *
 * This is an application framework for solving the problem of
 * predictive model selection of coupled models as presented in the 
 * following paper:
 * 
 * Gabriel Terejanu, Todd Oliver, Chris Simmons (2011). Application of 
 * Predictive Model Selection to Coupled Models. In Proceedings of the World 
 * Congress on Engineering and Computer Science 2011 Vol II, WCECS 2011, 
 * pp. 927-932.
 * 
 * The framework is built on top of statistical library QUESO 
 * (Quantification of Uncertainty for Estimation, Simulation and Optimization).
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *--------------------------------------------------------------------------
 *
 * forcing_SED_funcs_test.C
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <cmath>
#include <map>
#include <string>
#include <iostream>
#include "uqDefines.h"
#include "uqGslVector.h"
#include "uqVectorSpace.h"
#include "uqGslMatrix.h"
#include "testsCommon.h"
#include "forcingSimpleExponentialDecay.h"

uqFullEnvironmentClass* env;

//********************************************************************
// Definitions
//********************************************************************
bool test_SED_func();
bool test_SED_jac();
bool test_SED_qoi();

//********************************************************************
// main - SED tests
//********************************************************************
int main(int argc, char* argv[])
{

  MPI_Init(&argc,&argv);
  env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // tests SED_func
  TST_MACRO( test_SED_func(),
	     "simpleExponentialDecayForcing_func",
	     "" );

  // tests SED_jac
  TST_MACRO( test_SED_jac(),
	     "simpleExponentialDecayForcing_jac",
	     "" );

  // tests SED_qoi
  TST_MACRO( test_SED_qoi(),
	     "simpleExponentialDecayForcing_qoi",
	     "" );

  // Finalize environment
  delete env;
  MPI_Finalize();

  // no error has been found
  return 0;

}


//********************************************************************
// tests SED_func
//********************************************************************
bool test_SED_func()
{
  // inputs (scn & param)
  std::map<std::string, double> map_act_params;
  map_act_params[ __SED_PARAM_F0 ] = 2.0;
  map_act_params[ __SED_PARAM_TAU ] = 10.0;
  double scn_time = 0.0;
  
  // exact and compute value
  double exact_value = 2.0;
  double computed_value = simpleExponentialDecayForcing_func( scn_time,
							      (void *)&map_act_params );

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests SED_jac
//********************************************************************
bool test_SED_jac()
{
  // inputs (scn & param)
  std::map<std::string, double> map_act_params;
  map_act_params[ __SED_PARAM_F0 ] = 2.0;
  map_act_params[ __SED_PARAM_TAU ] = 1.0;
  double scn_time = 0.0;
  
  // exact and compute value
  double exact_value = -2.0;
  double computed_value = simpleExponentialDecayForcing_jac( scn_time,
							     (void *)&map_act_params );

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests SED_qoi
//********************************************************************
bool test_SED_qoi()
{

  // scenario
  double scn_time = __FORCING_INTEGRATION_PERIOD;

  // params
  double param_F0 = 2.0;
  double param_tau = 1.0;

  // map inputs (scn & param)

  // scenario does not matter here because we are using
  // __FORCING_INTEGRATION_PERIOD
  std::map<std::string, double> map_scn;
  map_scn[ __SED_SCN_TIME ] = scn_time;

  std::map<std::string, double> map_act_params;
  map_act_params[ __SED_PARAM_F0 ] = param_F0;
  map_act_params[ __SED_PARAM_TAU ] = param_tau;

  // prepare qoi vector
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>  
    qoiSpace( *env, "qoi", __SED_QOI_NUMBER, NULL );
      
  uqGslVectorClass qoiVec( qoiSpace.zeroVector() );
  
  // exact and compute value
  double exact_value = -param_F0*param_tau*
    (exp(-scn_time/param_tau)-1.0);

  simpleExponentialDecayForcing_qoi<uqGslVectorClass,uqGslMatrixClass>  
    ( qoiVec, map_scn, map_act_params );

  double computed_value = qoiVec[0];

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}


