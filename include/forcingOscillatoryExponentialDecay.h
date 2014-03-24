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
 * forcingOscillatoryExponentialDecay.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __FORCING_OSCILLATORY_EXPONENTIAL_DECAY_H
#define __FORCING_OSCILLATORY_EXPONENTIAL_DECAY_H

#include "forcingCommon.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


//********************************************************************
// oscillatoryExponentialDecayForcing : Definitions
//********************************************************************

// Observables
#define __OED_OBS_FORCE                   "obs_force"
#define __OED_OBS_FORCE_LN                true
// QoIs
#define __OED_QOI_NUMBER                  2
#define __OED_QOI_INTEGRATED_FORCE        "qoi_integrated_force"
#define __OED_QOI_MIN_DECAY_TIME          "qoi_min_decay_time"
// Scenarios
#define __OED_SCN_TIME                    "scn_time"
// Parameters
#define __OED_PARAM_F0                    "param_Force_OED_F0"
#define __OED_PARAM_TAU                   "param_Force_OED_tau"
#define __OED_PARAM_ALPHA                 "param_Force_OED_alpha"
#define __OED_PARAM_OMEGA                 "param_Force_OED_omega"


//********************************************************************
// oscillatoryExponentialDecayForcing_func
//********************************************************************
double oscillatoryExponentialDecayForcing_func( double x, 
						void* params )
{

  // inputs
  double scn_time = x;
  std::map<std::string, double> map_act_params = 
    *((std::map<std::string, double> *)params);

  // get parameters
  double param_F0 = map_act_params[ __OED_PARAM_F0 ];
  double param_tau = map_act_params[ __OED_PARAM_TAU ];
  double param_alpha = map_act_params[ __OED_PARAM_ALPHA ];
  double param_omega = map_act_params[ __OED_PARAM_OMEGA ];

  // outputs
  double f;

  // calculations
  f = param_F0 * exp( - scn_time / param_tau ) * 
    ( param_alpha * sin( param_omega * scn_time ) + 1.0 );
  
  return f;

}


//********************************************************************
// oscillatoryExponentialDecayForcing_jac
//********************************************************************
double oscillatoryExponentialDecayForcing_jac( double x, 
					       void* params )
{

  // inputs
  double scn_time = x;
  std::map<std::string, double> map_act_params = 
    *((std::map<std::string, double> *)params);

  // get parameters
  double param_F0 = map_act_params[ __OED_PARAM_F0 ];
  double param_tau = map_act_params[ __OED_PARAM_TAU ];
  double param_alpha = map_act_params[ __OED_PARAM_ALPHA ];
  double param_omega = map_act_params[ __OED_PARAM_OMEGA ];

  // outputs
  double dfdt;

  // calculations
  dfdt = -param_F0/param_tau * exp( - scn_time / param_tau ) * 
    ( param_alpha * sin( param_omega * scn_time ) + 1.0 ) +
    param_alpha*param_F0*param_omega*cos( param_omega * scn_time ) *
    exp( - scn_time / param_tau );
  
  return dfdt;

}


//********************************************************************
// oscillatoryExponentialDecayForcing_obs
//********************************************************************
template <class P_V, class P_M>
  void oscillatoryExponentialDecayForcing_obs( P_V& obsVec,
					       std::map<std::string, double> &map_act_params,
					       std::map<std::string, double> &map_data_point )
{

  // get scenario
  double scn_time = map_data_point[ __OED_SCN_TIME ];

  // compute the observable
  double f;
  f = oscillatoryExponentialDecayForcing_func( scn_time, &map_act_params );
  
  // save results
  obsVec[0] = f;

}


//********************************************************************
// oscillatoryExponentialDecayForcing_root
//********************************************************************
double oscillatoryExponentialDecayForcing_root( double x, 
					   void* params )
{

  double f = oscillatoryExponentialDecayForcing_func( x, params ) - 
    __FORCING_FORCE_DECAY_TOLERANCE;

  return f;

}

//********************************************************************
// oscillatoryExponentialDecayForcing_get_one_root
//********************************************************************
int oscillatoryExponentialDecayForcing_get_one_root
( std::map<std::string, double>& map_act_params,
  double time_low, double time_high, 
  double &ret_root)
{

  // bracket
  double x_lo = time_low;
  double x_hi = time_high;

  // gsl prep
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  gsl_function F;
  F.function = &oscillatoryExponentialDecayForcing_root;
  F.params = &map_act_params;
     
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc( T );
  gsl_root_fsolver_set( s, &F, x_lo, x_hi );
     
  // find root
  double r = 0;
  int status;
  int iter = 0, max_iter = 100;
  int ret_flag = GSL_FAILURE;

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate( s );

      r = gsl_root_fsolver_root( s );
      x_lo = gsl_root_fsolver_x_lower( s );
      x_hi = gsl_root_fsolver_x_upper( s );
      status = gsl_root_test_interval( x_lo, x_hi,
				       0, __FORCING_GSL_EPSROOT );

      /*
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
	      iter, x_lo, x_hi,
	      r, x_hi - x_lo);
      */

      if (status == GSL_SUCCESS) 
	{
	  ret_root = r;
	  ret_flag = GSL_SUCCESS;
	  break;
	}
     
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  // deallocate solver
  gsl_root_fsolver_free( s );

  return ret_flag;
}


//********************************************************************
// oscillatoryExponentialDecayForcing_qoi_get_time_decay
//********************************************************************
double oscillatoryExponentialDecayForcing_qoi_get_time_decay
  ( std::map<std::string, double>& map_act_params )
{

  double time_low = 0.0;
  double time_high = 0.0;
  double min_time_decay = __FORCING_INTEGRATION_PERIOD;

  double tmp_time_decay = 0.0;
  int err_code = GSL_SUCCESS;

  // get params
  double param_alpha = map_act_params[ __OED_PARAM_ALPHA ];
  double param_omega = map_act_params[ __OED_PARAM_OMEGA ];

  double qtr_period;
  int no_qtr_periods;
  if( ( param_alpha != 0.0 ) && ( param_omega != 0.0 ) )
    {
      qtr_period = PI / ( 2.0 * param_omega );
      no_qtr_periods = (int)(__FORCING_INTEGRATION_PERIOD / qtr_period);
    }
  else
    {
      qtr_period = __FORCING_INTEGRATION_PERIOD;
      no_qtr_periods = 1;
    }
  /*
  std::cout << "qtr_period = " << qtr_period << std::endl;
  std::cout << "no_qtr_periods = " << no_qtr_periods << std::endl;
  */

  // find the left most root
  for( unsigned int i = 0; i < no_qtr_periods; ++i )
    {

      // compute bracket
      time_low = qtr_period * (double)i;
      double f_high = 
	oscillatoryExponentialDecayForcing_root( time_low, (void *)&map_act_params );

      time_high = qtr_period * (double)(i+1);
      double f_low = 
	oscillatoryExponentialDecayForcing_root( time_high, (void *)&map_act_params );

      /*
      std::cout << "[ " << time_low << "(" << f_high << ") , " <<
	time_high << "(" << f_low << ") ]" << std::endl;
      */

      // check bracket
      if( ( f_high >= 0.0 ) && ( f_low <= 0.0 ) )
	{
	  err_code = oscillatoryExponentialDecayForcing_get_one_root( map_act_params,
								      time_low,
								      time_high,
								      tmp_time_decay );
	  if( err_code == GSL_SUCCESS )
	    {
	      min_time_decay = tmp_time_decay;
	      break;
	    }
	}

    }

    UQ_FATAL_TEST_MACRO
      ( min_time_decay == __FORCING_INTEGRATION_PERIOD,
	0,
	"oscillatoryExponentialDecayForcing : qoi_get_time_decay",
	"Could not find any roots" );

    //    std::cout << "min time decay = " << min_time_decay << std::endl;

    return min_time_decay;

}



//********************************************************************
// oscillatoryExponentialDecayForcing_qoi
//********************************************************************
template <class P_V, class P_M>
  void oscillatoryExponentialDecayForcing_qoi( P_V& qoiVec,
					       std::map<std::string, double>& map_scn,
					       std::map<std::string, double>& map_act_params )
{

  // get scenario (NOTE: here it does not matter, we use 
  // __FORCING_INTEGRATION_PERIOD for the entire integration period) 
  double scn_time = map_scn[ __OED_SCN_TIME ];

  // GSL prep
  gsl_integration_workspace* gsl_ws 
    = gsl_integration_workspace_alloc( __FORCING_QUADRATURE_INTERVALS );

  gsl_function gsl_F;
  gsl_F.function = &oscillatoryExponentialDecayForcing_func;
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
  qoiVec[1] = oscillatoryExponentialDecayForcing_qoi_get_time_decay( map_act_params );

}


#endif // __FORCING_OSCILLATORY_EXPONENTIAL_DECAY_H
