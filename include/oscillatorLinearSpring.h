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
 * oscillatorLinearSpring.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __OSCILLATOR_LINEAR_SPRING_H
#define __OSCILLATOR_LINEAR_SPRING_H

#include "oscillatorCommon.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "uqDefines.h"

//********************************************************************
// oscillatorLinearSpring : Definitions
//********************************************************************

// Observables
#define __OLS_OBS_KINETIC_ENERGY          "obs_kinetic_energy"
#define __OLS_OBS_KINETIC_ENERGY_LN       true
// QoIs
#define __OLS_QOI_NUMBER                  3
#define __OLS_QOI_TIME_TO_REST            "qoi_time_to_rest"
#define __OLS_QOI_MAX_DISPLACEMENT        "qoi_max_displacement"
#define __OLS_QOI_MAX_VELOCITY            "qoi_max_velocity"
// Scenarios
#define __OLS_SCN_TIME                    "scn_time"
// Parameters
#define __OLS_PARAM_C                     "param_Osc_OLS_c"
#define __OLS_PARAM_K10                   "param_Osc_OLS_k10"


//********************************************************************
// oscillatorLinearSpring_func
//********************************************************************
int oscillatorLinearSpring_func( double t, 
				 const double y[],
				 double f[],
				 void* params )
{

  // inputs
  std::map<std::string, double> map_act_params = 
    *((std::map<std::string, double> *)params);

  // get parameters
  double param_c = map_act_params[ __OLS_PARAM_C ];
  double param_k10 = map_act_params[ __OLS_PARAM_K10 ];

  // get force
  double force = getForce_func( t, map_act_params );

  // calculations
  f[0] = y[1];
  f[1] = force 
    - param_c/__OSCILLATOR_MASS*y[1] - param_k10/__OSCILLATOR_MASS*y[0];

  return GSL_SUCCESS;

}

//********************************************************************
// oscillatorLinearSpring_jac
//********************************************************************
int oscillatorLinearSpring_jac( double t, 
				const double y[],
				double *dfdy,
				double dfdt[],
				void* params )
{

  // inputs
  std::map<std::string, double> map_act_params = 
    *((std::map<std::string, double> *)params);

  // get parameters
  double param_c = map_act_params[ __OLS_PARAM_C ];
  double param_k10 = map_act_params[ __OLS_PARAM_K10 ];

  // gsl
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 2, 2);

  gsl_matrix * m = &dfdy_mat.matrix; 

  // get derivative force
  double deriv_force = getForce_jac( t, map_act_params );

  // calculations
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -param_k10/__OSCILLATOR_MASS );
  gsl_matrix_set (m, 1, 1, -param_c/__OSCILLATOR_MASS );
  
  dfdt[0] = 0.0;
  dfdt[1] = deriv_force;

  return GSL_SUCCESS;

}

//********************************************************************
// oscillatorLinearSpring_int
//********************************************************************
void oscillatorLinearSpring_int( double &final_y0, double &final_y1,
				 double &time_to_rest,
				 double &max_displacement,
				 double &max_velocity,
				 bool get_qoi,
				 std::map<std::string, double>& map_scn,
				 std::map<std::string, double>& map_act_params )
{

  // get scenario 
  double scn_time;
  if( get_qoi )
    scn_time = __OSCILLATOR_INTEGRATION_PERIOD;
  else
    scn_time = map_scn[ __OLS_SCN_TIME ];
  
  // get parameters
  double param_c = map_act_params[ __OLS_PARAM_C ];
  double param_k10 = map_act_params[ __OLS_PARAM_K10 ];

  // GSL prep
  unsigned int dim = 2;
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc( T, dim );
  gsl_odeiv_control * c = gsl_odeiv_control_y_new( __OSCILLATOR_GSL_EPSABS, 
						   __OSCILLATOR_GSL_EPSREL );
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc( dim );
  gsl_odeiv_system sys = { oscillatorLinearSpring_func, 
			   oscillatorLinearSpring_jac, 
			   dim, &map_act_params };
     
  // initial condition
  double y[2] = { __OSCILLATOR_IC_X1, __OSCILLATOR_IC_X2 };
  double t = 0.0, t1 = scn_time;
  double h = 1e-6;    // step-size
  
  // qois
  time_to_rest = t;
  bool time_to_rest_set = false;
  max_displacement = 0.0;
  max_velocity = 0.0;

  // GSL integration
  while (t < t1)
    {
      int status = gsl_odeiv_evolve_apply( e, c, s,
					   &sys, 
					   &t, t1,
					   &h, y );

      // printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
      UQ_FATAL_TEST_MACRO( status != GSL_SUCCESS,
			   0,
			   "oscillatorLinearSpring_int",
			   "The status of GSL integration != GSL_SUCCESS" );

      // QoI 1 : compute total energy = kinetic energy + potential energy
      double tmp_energy = pow( y[1], 2.0 )/2.0 + 
	param_k10/2.0 * pow( y[0], 2.0 );      

      if( ( tmp_energy < __OSCILLATOR_TOTAL_ENERGY_EPS ) && ( !time_to_rest_set ) )
	{
	  time_to_rest = t;
	  time_to_rest_set = true;
	}

      if( tmp_energy >= __OSCILLATOR_TOTAL_ENERGY_EPS ) 
	time_to_rest_set = false;

      // QoI 2 : maximum displacement
      if( fabs( y[0] ) > max_displacement )
	max_displacement = fabs( y[0] );

      // QoI 3 : maximum velocity
      if( fabs( y[1] ) > max_velocity )
	max_velocity = fabs( y[1] );

    }

  // save results
  final_y0 = y[0];
  final_y1 = y[1];

  // deallocate memory   
  gsl_odeiv_evolve_free( e );
  gsl_odeiv_control_free( c );
  gsl_odeiv_step_free( s );

}

//********************************************************************
// oscillatorLinearSpring_obs
//********************************************************************
template <class P_V, class P_M>
  void oscillatorLinearSpring_obs( P_V& obsVec,
				   std::map<std::string, double> &map_act_params,
				   std::map<std::string, double> &map_data_point )
{

  // outputs
  double final_y0;
  double final_y1;
  double time_to_rest;
  double max_displacement;
  double max_velocity;

  oscillatorLinearSpring_int( final_y0, final_y1,
			      time_to_rest,
			      max_displacement,
			      max_velocity,
			      false,
			      map_data_point,
			      map_act_params );

  // save results
  obsVec[0] = pow( final_y1, 2.0 ) / 2.0;                // KINETIC ENERGY

}

//********************************************************************
// oscillatorLinearSpring_qoi
//********************************************************************
template <class P_V, class P_M>
  void oscillatorLinearSpring_qoi( P_V& qoiVec,
				   std::map<std::string, double>& map_scn,
				   std::map<std::string, double>& map_act_params )
{

  // outputs
  double final_y0;
  double final_y1;
  double time_to_rest;
  double max_displacement;
  double max_velocity;

  oscillatorLinearSpring_int( final_y0, final_y1,
			      time_to_rest,
			      max_displacement,
			      max_velocity,
			      true,
			      map_scn,
			      map_act_params );

  // save results
  qoiVec[0] = time_to_rest;                // TIME TO REST
  qoiVec[1] = max_displacement;            // MAXIMUM DISPLACEMENT
  qoiVec[2] = max_velocity;                // MAXIMUM VELOCITY

}


#endif // __OSCILLATOR_LINEAR_SPRING_H
