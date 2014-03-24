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
 * modelOscillatorQuinticSpring.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __MODEL_OSCILLATOR_QUINTIC_SPRING_H
#define __MODEL_OSCILLATOR_QUINTIC_SPRING_H


#include "uqModel.h"
#include "uqDataSet.h"
#include "oscillatorQuinticSpring.h"

//********************************************************************
// modelOscillatorQuinticSpring Definitions
//********************************************************************

// Model Name
#define __OQS_NAME                        "OscOQS"
// Noise
#define __OQS_PARAM_SIGMA                 "param_Osc_OQS_sigma"
// #define __OQS_SIGMA                       0.05


//********************************************************************
// modelOscillatorQuinticSpring Class
//********************************************************************
template <class P_V, class P_M>
  class modelOscillatorQuinticSpring : public uqModel<P_V,P_M>
{

 public:

  modelOscillatorQuinticSpring( const uqBaseEnvironmentClass& env );
  modelOscillatorQuinticSpring( const uqBaseEnvironmentClass& env, 
			       std::string modelName );
  virtual ~modelOscillatorQuinticSpring();

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
// modelOscillatorQuinticSpring::modelOscillatorQuinticSpring
//********************************************************************
template <class P_V, class P_M>
  modelOscillatorQuinticSpring<P_V,P_M>::
  modelOscillatorQuinticSpring( const uqBaseEnvironmentClass& env )
  : uqModel<P_V,P_M>( env, __OQS_NAME )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// modelOscillatorQuinticSpring::modelOscillatorQuinticSpring
//********************************************************************
template <class P_V, class P_M>
  modelOscillatorQuinticSpring<P_V,P_M>::
  modelOscillatorQuinticSpring( const uqBaseEnvironmentClass& env, 
			       std::string modelName )
  : uqModel<P_V,P_M>( env, modelName )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// modelOscillatorQuinticSpring::~modelOscillatorQuinticSpring
//********************************************************************
template <class P_V, class P_M>
  modelOscillatorQuinticSpring<P_V,P_M>::~modelOscillatorQuinticSpring()
{}

//********************************************************************
// modelOscillatorQuinticSpring::my_forward_problem
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::
  my_forward_problem( P_V& predVec, 
		      std::map<std::string, double> &map_act_params, 
		      std::map<std::string, double> &map_data_point ) const
{

  oscillatorQuinticSpring_obs<P_V,P_M>( predVec, 
					map_act_params, 
					map_data_point );

}

//********************************************************************
// modelOscillatorQuinticSpring::my_noise_covariance_matrix
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::
  my_noise_covariance_matrix( P_M& covMat, 
			      const P_V& predVec,
			      std::map<std::string, double>& map_act_params ) const
{

  // 1. get it from a file

  // 2. compute it from the parameters
  double param_sigma = map_act_params[ __OQS_PARAM_SIGMA ];
  covMat(0,0) = pow( param_sigma, 2.0 ) + pow( __OSCILLATOR_MIN_SIGMA, 2.0 );

  // 3. compute it as a % from predicted values

  // 4. use a fixed value
  //  covMat(0,0) = pow( __OQS_SIGMA, 2.0 );

}

//********************************************************************
// modelOscillatorQuinticSpring::my_qoi_model
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::
  my_qoi_model( P_V& qoiPredVec,
		std::map<std::string, double>& map_scn,
		std::map<std::string, double>& map_act_params ) const
{

  oscillatorQuinticSpring_qoi<P_V,P_M>( qoiPredVec, 
					map_scn,
					map_act_params );

}

//********************************************************************
// modelOscillatorQuinticSpring::add_spec_params
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::add_spec_params_in_order()
{

  addParam( __OQS_PARAM_C );
  addParam( __OQS_PARAM_K50 );
  addParam( __OQS_PARAM_K52 );
  addParam( __OQS_PARAM_K54 );
  addParam( __OQS_PARAM_SIGMA );

}

//********************************************************************
// modelOscillatorQuinticSpring::add_spec_observables_in_order
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::add_spec_observables_in_order()
{

  addObservable( __OQS_OBS_KINETIC_ENERGY, __OQS_OBS_KINETIC_ENERGY_LN );

}

//********************************************************************
// modelOscillatorQuinticSpring::add_spec_qois_in_order
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::add_spec_qois_in_order()
{

  addQoI( __OQS_QOI_TIME_TO_REST );
  addQoI( __OQS_QOI_MAX_DISPLACEMENT );
  addQoI( __OQS_QOI_MAX_VELOCITY );

}

//********************************************************************
// modelOscillatorQuinticSpring::add_spec_scenarios
//********************************************************************
template <class P_V, class P_M>
  void modelOscillatorQuinticSpring<P_V,P_M>::add_spec_scenarios_in_order()
{

  addScn( __OQS_SCN_TIME );

}


#endif // __MODEL_OSCILLATOR_QUINTIC_SPRING_H

