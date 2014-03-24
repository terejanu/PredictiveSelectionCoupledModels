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
 * coupled_OLS_SED_funcs_test.C
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
#include "oscillatorLinearSpring.h"
#include "forcingSimpleExponentialDecay.h"
#include "modelForcingSimpleExponentialDecay.h"

uqFullEnvironmentClass* env;

//********************************************************************
// Definitions
//********************************************************************
bool test_OLS_SED_func();
bool test_OLS_SED_int_exact();
bool test_OLS_SED_int_energy();

//********************************************************************
// main - OLS - SED tests
//********************************************************************
int main(int argc, char* argv[])
{

  MPI_Init(&argc,&argv);
  env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // tests OLS_SED_func
  TST_MACRO( test_OLS_SED_func(),
	     "oscillatorLinearSpring_func",
	     "with forcing" );

  // tests OLS_SED_int_exact
  TST_MACRO( test_OLS_SED_int_exact(),
	     "oscillatorLinearSpring_int",
	     "exact solution with forcing" );

  // Finalize environment
  delete env;
  MPI_Finalize();

  // no error has been found
  return 0;

}

//********************************************************************
// tests OLS_SED_func
//********************************************************************
bool test_OLS_SED_func()
{

  // scenario
  double scn_time = 1.0;

  // params oscillator
  double param_c = 2.0;
  double param_k10 = param_c;

  // params forcing
  double param_F0 = 2.0;
  double param_tau = 1.0;

  // map inputs params (oscillator)
  std::map<std::string, double> map_act_params;
  map_act_params[ __OLS_PARAM_C ] = param_c;
  map_act_params[ __OLS_PARAM_K10 ] = param_k10;

  // forcing
  map_act_params[ __SED_PARAM_F0 ] = param_F0;
  map_act_params[ __SED_PARAM_TAU ] = param_tau;  

  // set coupling
  map_act_params[ __SED_NAME ] = 1.0;

  // outputs
  double ic_y[2] = {0.5, 1.5};
  double computed_f[2] = {0.0, 0.0};
  
  double exact_f[2] = 
    { ic_y[1], 
      param_F0*exp(-scn_time/param_tau) - param_c/__OSCILLATOR_MASS*(ic_y[0]+ic_y[1]) };

  // get values
  oscillatorLinearSpring_func( scn_time, 
			       ic_y,
			       computed_f,
			       (void*)&map_act_params );

  // compare the exact value with the computed one
    return ( ( std::fabs( exact_f[0] - computed_f[0] ) > DOUBLE_TOL ) ||
	     ( std::fabs( exact_f[1] - computed_f[1] ) > DOUBLE_TOL ) );

}


//********************************************************************
// tests OLS_SED_int_exact
//********************************************************************
bool test_OLS_SED_int_exact()
{

  // scenario
  double scn_time = 1.0;

  // params oscillator
  double param_c = 0.1;                   
  double param_k10 = 0.0;                 // to have a linear ODE in velocity

  // params forcing
  double param_F0 = 2.0;
  double param_tau = 1.0;

  // map inputs 
  std::map<std::string, double> map_scn;
  map_scn[ __OLS_SCN_TIME ] = scn_time;

  // map params oscillator
  std::map<std::string, double> map_act_params;
  map_act_params[ __OLS_PARAM_C ] = param_c;
  map_act_params[ __OLS_PARAM_K10 ] = param_k10;

  // forcing
  map_act_params[ __SED_PARAM_F0 ] = param_F0;
  map_act_params[ __SED_PARAM_TAU ] = param_tau;  

  // set coupling
  map_act_params[ __SED_NAME ] = 1.0;

  // outputs
  double final_y0;
  double final_y1;
  double time_to_rest;
  double max_displacement;
  double max_velocity;

  // get final velocity
  oscillatorLinearSpring_int( final_y0, final_y1,
			      time_to_rest,
			      max_displacement,
			      max_velocity,
			      false,
			      map_scn,
			      map_act_params );

  double my_const = __OSCILLATOR_IC_X2 + 
    param_tau * param_F0 / (__OSCILLATOR_MASS - param_tau * param_c );

  double pre_exp = param_tau * param_F0 / ( param_tau * param_c - __OSCILLATOR_MASS );

  double exact_y1 = exp( -param_c/__OSCILLATOR_MASS*scn_time ) *
    ( pre_exp * 
      exp( (param_tau*param_c-__OSCILLATOR_MASS)/(__OSCILLATOR_MASS*param_tau)*scn_time  ) 
      + my_const );

  //  std::cout << "final_y1 = " << final_y1 << std::endl;
  //  std::cout << "exact_y1 = " << exact_y1 << std::endl;

  // compare the exact value with the computed one
  return ( std::fabs( exact_y1 - final_y1 ) > DOUBLE_TOL ) ;

}

