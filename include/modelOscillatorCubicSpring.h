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
 * modelOscillatorCubicSpring.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __MODEL_OSCILLATOR_CUBIC_SPRING_H
#define __MODEL_OSCILLATOR_CUBIC_SPRING_H


#include "uqModel.h"
#include "uqDataSet.h"
#include "oscillatorCubicSpring.h"


//********************************************************************
// modelOscillatorCubicSpring Definitions
//********************************************************************

// Model name
#define __OCS_NAME                        "OscOCS"
// Noise
#define __OCS_PARAM_SIGMA                 "param_Osc_OCS_sigma"
// #define __OCS_SIGMA                       0.05


//********************************************************************
// modelOscillatorCubicSpring Class
//********************************************************************
template <class P_V, class P_M>
  class modelOscillatorCubicSpring : public uqModel<P_V,P_M>
{

 public:

  modelOscillatorCubicSpring( const uqBaseEnvironmentClass& env );
  modelOscillatorCubicSpring( const uqBaseEnvironmentClass& env, 
			       std::string modelName );
  virtual ~modelOscillatorCubicSpring();

  using uqModel<P_V,P_M>::getActualParameters;
  using uqModel<P_V,P_M>::getNumberObservables;
  
 private:

  using uqModel<P_V,P_M>::m_prefix;
  using uqModel<P_V,P_M>::m_env;
  using uqModel<P_V,P_M>::m_num_params;
  using uqModel<P_V,P_M>::m_num_params_set;
  using uqModel<P_V,P_M>::m_set_params;

  using uqModel<P_V,P_M>::addParam;
  using uqModel<P_V,P_M>::addObservable;
  using uqModel<P_V,P_M>::addQoI;
  using uqModel<P_V,P_M>::addScn;
  using uqModel<P_V,P_M>::get_act_val_param;
  using uqModel<P_V,P_M>::compute_obs_pred_diff;
  using uqModel<P_V,P_M>::get_correction_term;

  void add_spec_params_in_order();
  void add_spec_observables_in_order();
  void add_spec_qois_in_order();
  void add_spec_scenarios_in_order();

  void my_forward_problem( P_V& predVec, 
			   std::map<std::string, double>& map_act_params, 
			   std::map<std::string, double>& map_data_point ) const;
  void my_noise_covariance_matrix( P_M& covMat, 
				   const P_V& predVec,
				   std::map<std::string, double>& map_act_params ) const;
  void my_qoi_model( P_V& qoiPredVec,
		     std::map<std::string, double>& map_scn,
		     std::map<std::string, double>& map_act_params ) const;

};

//********************************************************************
// modelOscillatorCubicSpring::modelOscillatorCubicSpring
//********************************************************************
template <class P_V, class P_M>
  modelOscillatorCubicSpring<P_V,P_M>::
  modelOscillatorCubicSpring( const uqBaseEnvironmentClass& env )
  : uqModel<P_V,P_M>( env, __OCS_NAME )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// modelOscillatorCubicSpring::modelOscillatorCubicSpring
//********************************************************************
template <class P_V, class P_M>
  modelOscillatorCubicSpring<P_V,P_M>::
  modelOscillatorCubicSpring( const uqBaseEnvironmentClass& env, 
			       std::string modelName )
  : uqModel<P_V,P_M>( env, modelName )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// modelOscillatorCubicSpring::~modelOscillatorCubicSpring
//********************************************************************
template <class P_V, class P_M>
  modelOscillatorCubicSpring<P_V,P_M>::~modelOscillatorCubicSpring()
{}

//********************************************************************
// modelOscillatorCubicSpring::my_forward_problem
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::
  my_forward_problem( P_V& predVec, 
		      std::map<std::string, double> &map_act_params, 
		      std::map<std::string, double> &map_data_point ) const
{

  oscillatorCubicSpring_obs<P_V,P_M>( predVec, 
				      map_act_params, 
				      map_data_point );

}

//********************************************************************
// modelOscillatorCubicSpring::my_noise_covariance_matrix
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::
  my_noise_covariance_matrix( P_M& covMat, 
			      const P_V& predVec,
			      std::map<std::string, double>& map_act_params ) const
{

  // 1. get it from a file

  // 2. compute it from the parameters
  double param_sigma = map_act_params[ __OCS_PARAM_SIGMA ];
  covMat(0,0) = pow( param_sigma, 2.0 ) + pow( __OSCILLATOR_MIN_SIGMA, 2.0 );

  // 3. compute it as a % from predicted values

  // 4. use a fixed value
  //  covMat(0,0) = pow( __OCS_SIGMA, 2.0 );

}

//********************************************************************
// modelOscillatorCubicSpring::my_qoi_model
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::
  my_qoi_model( P_V& qoiPredVec,
		std::map<std::string, double>& map_scn,
		std::map<std::string, double>& map_act_params ) const
{

  oscillatorCubicSpring_qoi<P_V,P_M>( qoiPredVec, 
				      map_scn,
				      map_act_params );

}

//********************************************************************
// modelOscillatorCubicSpring::add_spec_params
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::add_spec_params_in_order()
{

  addParam( __OCS_PARAM_C );
  addParam( __OCS_PARAM_K30 );
  addParam( __OCS_PARAM_K32 );
  addParam( __OCS_PARAM_SIGMA );

}

//********************************************************************
// modelOscillatorCubicSpring::add_spec_observables_in_order
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::add_spec_observables_in_order()
{

  addObservable( __OCS_OBS_KINETIC_ENERGY, __OCS_OBS_KINETIC_ENERGY_LN );

}

//********************************************************************
// modelOscillatorCubicSpring::add_spec_qois_in_order
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::add_spec_qois_in_order()
{

  addQoI( __OCS_QOI_TIME_TO_REST );
  addQoI( __OCS_QOI_MAX_DISPLACEMENT );
  addQoI( __OCS_QOI_MAX_VELOCITY );

}

//********************************************************************
// modelOscillatorCubicSpring::add_spec_scenarios
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorCubicSpring<P_V,P_M>::add_spec_scenarios_in_order()
{

  addScn( __OCS_SCN_TIME );

}


#endif // __MODEL_OSCILLATOR_CUBIC_SPRING_H

