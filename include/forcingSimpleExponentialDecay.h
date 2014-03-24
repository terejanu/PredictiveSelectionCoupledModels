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
 * forcingSimpleExponentialDecay.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __FORCING_SIMPLE_EXPONENTIAL_DECAY_H
#define __FORCING_SIMPLE_EXPONENTIAL_DECAY_H

#include "forcingCommon.h"
#include <gsl/gsl_integration.h>

//********************************************************************
// simpleExponentialDecayForcing : Definitions
//********************************************************************

// Observables
#define __SED_OBS_FORCE                   "obs_force"
#define __SED_OBS_FORCE_LN                true
// QoIs
#define __SED_QOI_NUMBER                  2
#define __SED_QOI_INTEGRATED_FORCE        "qoi_integrated_force"
#define __SED_QOI_MIN_DECAY_TIME          "qoi_min_decay_time"
// Scenarios
#define __SED_SCN_TIME                    "scn_time"
// Parameters
#define __SED_PARAM_F0                    "param_Force_SED_F0"
#define __SED_PARAM_TAU                   "param_Force_SED_tau"


//********************************************************************
// simpleExponentialDecayForcing_func
//********************************************************************
double simpleExponentialDecayForcing_func( double x, 
					   void* params )
{

  // inputs
  double scn_time = x;
  std::map<std::string, double> map_act_params = 
    *((std::map<std::string, double> *)params);

  // get parameters
  double param_F0 = map_act_params[ __SED_PARAM_F0 ];
  double param_tau = map_act_params[ __SED_PARAM_TAU ];

  // outputs
  double f;

  // calculations
  f = param_F0 * exp( - scn_time / param_tau );

  return f;

}

//********************************************************************
// simpleExponentialDecayForcing_jac
//********************************************************************
double simpleExponentialDecayForcing_jac( double x, 
					  void* params )
{

  // inputs
  double scn_time = x;
  std::map<std::string, double> map_act_params = 
    *((std::map<std::string, double> *)params);

  // get parameters
  double param_F0 = map_act_params[ __SED_PARAM_F0 ];
  double param_tau = map_act_params[ __SED_PARAM_TAU ];

  // outputs
  double dfdt;

  // calculations
  dfdt = - param_F0 / param_tau * exp( - scn_time / param_tau );

  return dfdt;

}

//********************************************************************
// simpleExponentialDecayForcing_obs
//********************************************************************
template <class P_V, class P_M>
  void simpleExponentialDecayForcing_obs( P_V& obsVec,
					  std::map<std::string, double> &map_act_params,
					  std::map<std::string, double> &map_data_point )
{

  // get scenario
  double scn_time = map_data_point[ __SED_SCN_TIME ];

  // compute the observable
  double f;
  f = simpleExponentialDecayForcing_func( scn_time, &map_act_params );
  
  // save results
  obsVec[0] = f;

}


//********************************************************************
// simpleExponentialDecayForcing_qoi
//********************************************************************
template <class P_V, class P_M>
  void simpleExponentialDecayForcing_qoi( P_V& qoiVec,
					  std::map<std::string, double>& map_scn,
					  std::map<std::string, double>& map_act_params )
{

  // get scenario (NOTE: here it does not matter, we use 
  // __FORCING_INTEGRATION_PERIOD for the entire integration period) 
  double scn_time = map_scn[ __SED_SCN_TIME ];


  // GSL prep
  gsl_integration_workspace* gsl_ws 
    = gsl_integration_workspace_alloc( __FORCING_QUADRATURE_INTERVALS );

  gsl_function gsl_F;
  gsl_F.function = &simpleExponentialDecayForcing_func;
  gsl_F.params = &map_act_params;
  double gsl_result, gsl_error;

  // GSL integration
  gsl_integration_qags( &gsl_F, 
			0.0, __FORCING_INTEGRATION_PERIOD, 
			__FORCING_GSL_EPSABS, __FORCING_GSL_EPSREL, 
			__FORCING_QUADRATURE_INTERVALS,
			gsl_ws, &gsl_result, &gsl_error );

  // save results ( qoi: integrated force )
  qoiVec[0] = gsl_result;

  // deallocate memory
  gsl_integration_workspace_free( gsl_ws );

  // save results ( qoi: min decay time )
  double param_F0 = map_act_params[ __SED_PARAM_F0 ];
  double param_tau = map_act_params[ __SED_PARAM_TAU ];

  qoiVec[1] = - param_tau * log( __FORCING_FORCE_DECAY_TOLERANCE / param_F0 );


}


#endif // __FORCING_SIMPLE_EXPONENTIAL_DECAY_H
