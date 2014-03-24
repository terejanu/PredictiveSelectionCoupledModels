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
 * forcing_OLD_funcs_test.C
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
#include "forcingOscillatoryLinearDecay.h"

uqFullEnvironmentClass* env;

//********************************************************************
// Definitions
//********************************************************************
bool test_OLD_func_t_le_tau();
bool test_OLD_func_t_g_tau();
bool test_OLD_jac_t_le_tau();
bool test_OLD_jac_t_g_tau();
bool test_OLD_qoi_t_le_tau();
bool test_OLD_qoi_t_g_tau();

//********************************************************************
// main - OLD tests
//********************************************************************
int main(int argc, char* argv[])
{

  MPI_Init(&argc,&argv);
  env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // tests OLD_func ( t <= tau )
  TST_MACRO( test_OLD_func_t_le_tau(),
	     "oscillatoryLinearDecayForcing_func",
	     "t <= tau" );

  // tests OLD_func ( t > tau )
  TST_MACRO( test_OLD_func_t_g_tau(),
	     "oscillatoryLinearDecayForcing_func",
	     "t > tau" );

  // tests OLD_jac ( t <= tau )
  TST_MACRO( test_OLD_jac_t_le_tau(),
	     "oscillatoryLinearDecayForcing_jac",
	     "t <= tau" );

  // tests OLD_jac ( t > tau )
  TST_MACRO( test_OLD_jac_t_g_tau(),
	     "oscillatoryLinearDecayForcing_jac",
	     "t > tau" );

  std::cout << "enter test 1" << std::endl;
  // tests OLD_qoi ( t <= tau )
  TST_MACRO( test_OLD_qoi_t_le_tau(),
	     "oscillatoryLinearDecayForcing_qoi",
	     "t <= tau" );

  std::cout << "enter test 2" << std::endl;
  // tests OLD_qoi ( t > tau )
  TST_MACRO( test_OLD_qoi_t_g_tau(),
	     "oscillatoryLinearDecayForcing_qoi",
	     "t > tau" );

  // Finalize environment
  delete env;
  MPI_Finalize();

  // no error has been found
  return 0;

}


//********************************************************************
// tests OLD_func ( t <= tau )
//********************************************************************
bool test_OLD_func_t_le_tau()
{
  // inputs (scn & param)
  std::map<std::string, double> map_act_params;
  map_act_params[ __OLD_PARAM_F0 ] = 2.0;
  map_act_params[ __OLD_PARAM_TAU ] = 10.0;
  map_act_params[ __OLD_PARAM_ALPHA ] = 1.0;
  map_act_params[ __OLD_PARAM_OMEGA ] = 0.0;
  double scn_time = 5.0;
  
  // exact and compute value
  double exact_value = 1.0;
  double computed_value = oscillatoryLinearDecayForcing_func( scn_time,
							      (void *)&map_act_params );

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests OLD_func ( t > tau )
//********************************************************************
bool test_OLD_func_t_g_tau()
{
  // inputs (scn & param)
  std::map<std::string, double> map_act_params;
  map_act_params[ __OLD_PARAM_F0 ] = 2.0;
  map_act_params[ __OLD_PARAM_TAU ] = 10.0;
  map_act_params[ __OLD_PARAM_ALPHA ] = 1.0;
  map_act_params[ __OLD_PARAM_OMEGA ] = 0.0;
  double scn_time = 15.0;
  
  // exact and compute value
  double exact_value = 0.0;
  double computed_value = oscillatoryLinearDecayForcing_func( scn_time,
							      (void *)&map_act_params );

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests OLD_jac ( t <= tau )
//********************************************************************
bool test_OLD_jac_t_le_tau()
{
  // inputs (scn & param)
  std::map<std::string, double> map_act_params;
  map_act_params[ __OLD_PARAM_F0 ] = 2.0;
  map_act_params[ __OLD_PARAM_TAU ] = 1.0;
  map_act_params[ __OLD_PARAM_ALPHA ] = 1.0;
  map_act_params[ __OLD_PARAM_OMEGA ] = 0.0;
  double scn_time = 0.5;
  
  // exact and compute value
  double exact_value = -2.0;
  double computed_value = oscillatoryLinearDecayForcing_jac( scn_time,
							     (void *)&map_act_params );

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests OLD_jac ( t > tau )
//********************************************************************
bool test_OLD_jac_t_g_tau()
{
  // inputs (scn & param)
  std::map<std::string, double> map_act_params;
  map_act_params[ __OLD_PARAM_F0 ] = 2.0;
  map_act_params[ __OLD_PARAM_TAU ] = 10.0;
  map_act_params[ __OLD_PARAM_ALPHA ] = 1.0;
  map_act_params[ __OLD_PARAM_OMEGA ] = 0.0;
  double scn_time = 15.0;
  
  // exact and compute value
  double exact_value = 0.0;
  double computed_value = oscillatoryLinearDecayForcing_jac( scn_time,
							     (void *)&map_act_params );

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests OLD_qoi ( t <= tau )
//********************************************************************
bool test_OLD_qoi_t_le_tau()
{

  // scenario
  double scn_time = __FORCING_INTEGRATION_PERIOD;

  // params
  double param_F0 = 0.8;
  double param_tau = 2.0 * scn_time;
  double param_alpha = 1.0;
  double param_omega = 0.0;

  // map inputs (scn & param)

  // scenario does not matter here because we are using
  // __FORCING_INTEGRATION_PERIOD
  std::map<std::string, double> map_scn;
  map_scn[ __OLD_SCN_TIME ] = scn_time;

  std::map<std::string, double> map_act_params;
  map_act_params[ __OLD_PARAM_F0 ] = param_F0;
  map_act_params[ __OLD_PARAM_TAU ] = param_tau;
  map_act_params[ __OLD_PARAM_ALPHA ] = param_alpha;
  map_act_params[ __OLD_PARAM_OMEGA ] = param_omega;

  // prepare qoi vector
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>  
    qoiSpace( *env, "qoi", __OLD_QOI_NUMBER, NULL );
      
  uqGslVectorClass qoiVec( qoiSpace.zeroVector() );
  
  // exact and compute value
  double exact_value = param_F0*(scn_time - pow( scn_time, 2.0)/(2.0*param_tau));

  oscillatoryLinearDecayForcing_qoi<uqGslVectorClass,uqGslMatrixClass>  
    ( qoiVec, map_scn, map_act_params );

  double computed_value = qoiVec[0];

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}

//********************************************************************
// tests OLD_qoi ( t > tau )
//********************************************************************
bool test_OLD_qoi_t_g_tau()
{

  // scenario
  double scn_time = __FORCING_INTEGRATION_PERIOD;

  // params
  double param_F0 = 2.0;
  double param_tau = scn_time/2.0;
  double param_alpha = 1.0;
  double param_omega = 0.0;

  // map inputs (scn & param)

  // scenario does not matter here because we are using
  // __FORCING_INTEGRATION_PERIOD
  std::map<std::string, double> map_scn;
  map_scn[ __OLD_SCN_TIME ] = scn_time;

  std::map<std::string, double> map_act_params;
  map_act_params[ __OLD_PARAM_F0 ] = param_F0;
  map_act_params[ __OLD_PARAM_TAU ] = param_tau;
  map_act_params[ __OLD_PARAM_ALPHA ] = param_alpha;
  map_act_params[ __OLD_PARAM_OMEGA ] = param_omega;

  // prepare qoi vector
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>  
    qoiSpace( *env, "qoi", __OLD_QOI_NUMBER, NULL );
      
  uqGslVectorClass qoiVec( qoiSpace.zeroVector() );
  
  // exact and compute value
  double exact_value = param_F0*param_tau/2.0;

  oscillatoryLinearDecayForcing_qoi<uqGslVectorClass,uqGslMatrixClass>  
    ( qoiVec, map_scn, map_act_params );

  double computed_value = qoiVec[0];

  // compare the exact value with the computed one
  return ( std::fabs( exact_value - computed_value ) > DOUBLE_TOL ) ;
}
