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
 * oscillator_OQS_funcs_test.C
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
#include "oscillatorQuinticSpring.h"

uqFullEnvironmentClass* env;

//********************************************************************
// Definitions
//********************************************************************
bool test_OQS_func();
bool test_OQS_int_exact();
bool test_OQS_int_energy();

//********************************************************************
// main - OQS tests
//********************************************************************
int main(int argc, char* argv[])
{

  MPI_Init(&argc,&argv);
  env = new uqFullEnvironmentClass(MPI_COMM_WORLD,argv[1],"",NULL);

  // tests OQS_func
  TST_MACRO( test_OQS_func(),
	     "oscillatorQuinticSpring_func",
	     "" );

  // tests OQS_int_exact
  TST_MACRO( test_OQS_int_exact(),
	     "oscillatorQuinticSpring_int",
	     "exact solution" );

  // tests OQS_int_energy
  TST_MACRO( test_OQS_int_energy(),
	     "oscillatorQuinticSpring_int",
	     "const energy" );

  // Finalize environment
  delete env;
  MPI_Finalize();

  // no error has been found
  return 0;

}

//********************************************************************
// tests OQS_func
//********************************************************************
bool test_OQS_func()
{

  // scenario
  double scn_time = 1.0;

  // params
  double param_c = 2.0;
  double param_k50 = param_c;
  double param_k52 = param_c;
  double param_k54 = param_c;

  // map inputs params
  std::map<std::string, double> map_act_params;
  map_act_params[ __OQS_PARAM_C ] = param_c;
  map_act_params[ __OQS_PARAM_K50 ] = param_k50;
  map_act_params[ __OQS_PARAM_K52 ] = param_k52;
  map_act_params[ __OQS_PARAM_K54 ] = param_k54;

  // outputs
  double ic_y[2] = {0.5, 1.5};
  double computed_f[2] = {0.0, 0.0};
  
  double exact_f[2] = { ic_y[1], -param_c/__OSCILLATOR_MASS*(ic_y[0]+ic_y[1]+pow(ic_y[0],3.0)+pow(ic_y[0],5.0)) };

  // get values
  oscillatorQuinticSpring_func( scn_time, 
				ic_y,
				computed_f,
				(void*)&map_act_params );

  // compare the exact value with the computed one
  return ( ( std::fabs( exact_f[0] - computed_f[0] ) > DOUBLE_TOL ) ||
	   ( std::fabs( exact_f[1] - computed_f[1] ) > DOUBLE_TOL ) );

}


//********************************************************************
// tests OQS_int_exact
//********************************************************************
bool test_OQS_int_exact()
{

  // scenario
  double scn_time = 1.0;

  // params
  double param_c = 0.1;                   
  double param_k50 = 0.0;                 // to have a linear ODE in velocity
  double param_k52 = 0.0;                 // to have a linear ODE in velocity
  double param_k54 = 0.0;                 // to have a linear ODE in velocity

  // map inputs (scn & param)
  std::map<std::string, double> map_scn;
  map_scn[ __OQS_SCN_TIME ] = scn_time;

  std::map<std::string, double> map_act_params;
  map_act_params[ __OQS_PARAM_C ] = param_c;
  map_act_params[ __OQS_PARAM_K50 ] = param_k50;
  map_act_params[ __OQS_PARAM_K52 ] = param_k52;
  map_act_params[ __OQS_PARAM_K54 ] = param_k54;

  // outputs
  double final_y0;
  double final_y1;
  double time_to_rest;
  double max_displacement;
  double max_velocity;

  // get final velocity
  oscillatorQuinticSpring_int( final_y0, final_y1,
			       time_to_rest,
			       max_displacement,
			       max_velocity,
			       false,
			       map_scn,
			       map_act_params );

  double exact_y1 = __OSCILLATOR_IC_X2 * exp( -param_c/__OSCILLATOR_MASS*scn_time );

  // compare the exact value with the computed one
  return ( std::fabs( exact_y1 - final_y1 ) > DOUBLE_TOL ) ;

}


//********************************************************************
// tests OQS_int_energy
// NOTE: when we have no damping the energy has to be constant
//********************************************************************
bool test_OQS_int_energy()
{

  // params
  double param_c = 0.0;                   // no damping
  double param_k50 = 4.0;
  double param_k52 = -5.0;
  double param_k54 = 1.0;  

  // map inputs (scn & param)
  std::map<std::string, double> map_scn;
  std::map<std::string, double> map_act_params;
  map_act_params[ __OQS_PARAM_C ] = param_c;
  map_act_params[ __OQS_PARAM_K50 ] = param_k50;
  map_act_params[ __OQS_PARAM_K52 ] = param_k52;
  map_act_params[ __OQS_PARAM_K54 ] = param_k54;

  // outputs
  double final_y0;
  double final_y1;
  double time_to_rest;
  double max_displacement;
  double max_velocity;

  // get the energy at time: sc_time_1
  double scn_time_1 = 1.0;
  map_scn[ __OQS_SCN_TIME ] = scn_time_1;

  oscillatorQuinticSpring_int( final_y0, final_y1,
			       time_to_rest,
			       max_displacement,
			       max_velocity,
			       false,
			       map_scn,
			       map_act_params );

  double energy1 = pow( final_y1, 2.0 ) / 2.0 +
    pow( final_y0, 2.0 ) / 2.0 * param_k50 +
    pow( final_y0, 4.0 ) / 4.0 * param_k52 +
    pow( final_y0, 6.0 ) / 6.0 * param_k54;

  // get the energy at time: sc_time_2
  double scn_time_2 = 1.05;
  map_scn[ __OQS_SCN_TIME ] = scn_time_2;

  oscillatorQuinticSpring_int( final_y0, final_y1,
			       time_to_rest,
			       max_displacement,
			       max_velocity,
			       false,
			       map_scn,
			       map_act_params );

  double energy2 = pow( final_y1, 2.0 ) / 2.0 +
    pow( final_y0, 2.0 ) / 2.0 * param_k50 +
    pow( final_y0, 4.0 ) / 4.0 * param_k52 +
    pow( final_y0, 6.0 ) / 6.0 * param_k54;

  //  std::cout << "E1 = " << energy1 << std::endl;
  //  std::cout << "E2 = " << energy2 << std::endl;

  // compare the exact value with the computed one
  return ( std::fabs( energy1 - energy2 ) > DOUBLE_TOL ) ;

}
