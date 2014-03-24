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
 * oscillatorCommon.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __OSCILLATOR_COMMON_H
#define __OSCILLATOR_COMMON_H


#include "uqDefines.h"
#include "modelForcingSimpleLinearDecay.h"
#include "forcingSimpleLinearDecay.h"
#include "modelForcingSimpleExponentialDecay.h"
#include "forcingSimpleExponentialDecay.h"
#include "modelForcingOscillatoryLinearDecay.h"
#include "forcingOscillatoryLinearDecay.h"
#include "modelForcingOscillatoryExponentialDecay.h"
#include "forcingOscillatoryExponentialDecay.h"


#define __OSCILLATOR_GSL_EPSABS                   1e-8
#define __OSCILLATOR_GSL_EPSREL                   1e-8

#define __OSCILLATOR_MASS                         1.0          // kg
#define __OSCILLATOR_IC_X1                        0.5
#define __OSCILLATOR_IC_X2                        0.5

#define __OSCILLATOR_MIN_SIGMA                    0.1          // known instrument error

#define __OSCILLATOR_TOTAL_ENERGY_EPS             0.1
#define __OSCILLATOR_INTEGRATION_PERIOD           100.0     

//********************************************************************
// getForce_func
//********************************************************************
double getForce_func( double t,
		      std::map<std::string, double> &map_act_params )
{

  std::map<std::string, double>::const_iterator it;

  // ForceSLD
  it = map_act_params.find( __SLD_NAME );
  if( it != map_act_params.end() )
    {
      return simpleLinearDecayForcing_func( t, &map_act_params );
    }

  // ForceSED
  it = map_act_params.find( __SED_NAME );
  if( it != map_act_params.end() )
    {
      return simpleExponentialDecayForcing_func( t, &map_act_params );
    }

  // ForceOLD
  it = map_act_params.find( __OLD_NAME );
  if( it != map_act_params.end() )
    {
      return oscillatoryLinearDecayForcing_func( t, &map_act_params );
    }

  // ForceOED
  it = map_act_params.find( __OED_NAME );
  if( it != map_act_params.end() )
    {
      return oscillatoryExponentialDecayForcing_func( t, &map_act_params );
    }

  // key is not assigned => no forcing
  return 0.0;

}


//********************************************************************
// getForce_jac
//********************************************************************
double getForce_jac( double t,
		     std::map<std::string, double> &map_act_params )
{

  std::map<std::string, double>::const_iterator it;

  // ForceSLD
  it = map_act_params.find( __SLD_NAME );
  if( it != map_act_params.end() )
    {
      return simpleLinearDecayForcing_jac( t, &map_act_params );
    }

  // ForceSED
  it = map_act_params.find( __SED_NAME );
  if( it != map_act_params.end() )
    {
      return simpleExponentialDecayForcing_jac( t, &map_act_params );
    }

  // ForceOLD
  it = map_act_params.find( __OLD_NAME );
  if( it != map_act_params.end() )
    {
      return oscillatoryLinearDecayForcing_jac( t, &map_act_params );
    }

  // ForceOED
  it = map_act_params.find( __OED_NAME );
  if( it != map_act_params.end() )
    {
      return oscillatoryExponentialDecayForcing_jac( t, &map_act_params );
    }

  // key is not assigned => no forcing
  return 0.0;

}


#endif // __OSCILLATOR_COMMON_H
