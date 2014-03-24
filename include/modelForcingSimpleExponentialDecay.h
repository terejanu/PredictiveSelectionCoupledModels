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
 * modelForcingSimpleExponentialDecay.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __MODEL_FORCING_SIMPLE_EXPONENTIAL_DECAY_H
#define __MODEL_FORCING_SIMPLE_EXPONENTIAL_DECAY_H


#include "uqModel.h"
#include "uqDataSet.h"
#include "forcingSimpleExponentialDecay.h"

//********************************************************************
// modelForcingSimpleExponentialDecay Definitions
//********************************************************************

// Model name
#define __SED_NAME                        "ForceSED"
// Noise
#define __SED_PARAM_SIGMA                 "param_Force_SED_sigma"
// #define __SED_SIGMA                       0.05


//********************************************************************
// modelForcingSimpleExponentialDecay Class
//********************************************************************
template <class P_V, class P_M>
  class modelForcingSimpleExponentialDecay : public uqModel<P_V,P_M>
{

 public:

  modelForcingSimpleExponentialDecay( const uqBaseEnvironmentClass& env );
  modelForcingSimpleExponentialDecay( const uqBaseEnvironmentClass& env, 
				      std::string modelName );
  virtual ~modelForcingSimpleExponentialDecay();

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
// modelForcingSimpleExponentialDecay::modelForcingSimpleExponentialDecay
//********************************************************************
template <class P_V, class P_M>
  modelForcingSimpleExponentialDecay<P_V,P_M>::
  modelForcingSimpleExponentialDecay( const uqBaseEnvironmentClass& env )
  : uqModel<P_V,P_M>( env, __SED_NAME )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// modelForcingSimpleExponentialDecay::modelForcingSimpleExponentialDecay
//********************************************************************
template <class P_V, class P_M>
  modelForcingSimpleExponentialDecay<P_V,P_M>::
  modelForcingSimpleExponentialDecay( const uqBaseEnvironmentClass& env, 
				      std::string modelName )
  : uqModel<P_V,P_M>( env, modelName )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// modelForcingSimpleExponentialDecay::~modelForcingSimpleExponentialDecay
//********************************************************************
template <class P_V, class P_M>
  modelForcingSimpleExponentialDecay<P_V,P_M>::~modelForcingSimpleExponentialDecay()
{}


//********************************************************************
// modelForcingSimpleExponentialDecay::my_forward_problem
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::
  my_forward_problem( P_V& predVec, 
		      std::map<std::string, double> &map_act_params, 
		      std::map<std::string, double> &map_data_point ) const
{

  simpleExponentialDecayForcing_obs<P_V,P_M>( predVec, 
					      map_act_params,
					      map_data_point );

}


//********************************************************************
// modelForcingSimpleExponentialDecay::my_noise_covariance_matrix
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::
  my_noise_covariance_matrix( P_M& covMat, 
			      const P_V& predVec,
			      std::map<std::string, double>& map_act_params ) const
{

  // 1. get it from a file

  // 2. compute it from the parameters
  double param_sigma = map_act_params[ __SED_PARAM_SIGMA ];
  covMat(0,0) = pow( param_sigma, 2.0 ) + pow( __FORCING_MIN_SIGMA, 2.0 );

  // 3. compute it as a % from predicted values
  //  double param_sigma = map_act_params[ __SED_PARAM_SIGMA ];
  //  param_sigma = param_sigma * predVec[0];
  //  covMat(0,0) = pow( param_sigma, 2.0 );

  // 4. use a fixed value
  //  covMat(0,0) = pow( __SED_SIGMA, 2.0 );

}

//********************************************************************
// modelForcingSimpleExponentialDecay::my_qoi_model
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::
  my_qoi_model( P_V& qoiPredVec,
		std::map<std::string, double>& map_scn,
		std::map<std::string, double>& map_act_params ) const
{

  simpleExponentialDecayForcing_qoi<P_V,P_M>( qoiPredVec, 
					      map_scn,
					      map_act_params );

}


//********************************************************************
// modelForcingSimpleExponentialDecay::add_spec_params
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::add_spec_params_in_order()
{

  addParam( __SED_PARAM_F0 );
  addParam( __SED_PARAM_TAU );
  addParam( __SED_PARAM_SIGMA );

}

//********************************************************************
// modelForcingSimpleExponentialDecay::add_spec_observables_in_order
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::add_spec_observables_in_order()
{

  addObservable( __SED_OBS_FORCE, __SED_OBS_FORCE_LN );

}

//********************************************************************
// modelForcingSimpleExponentialDecay::add_spec_qois_in_order
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::add_spec_qois_in_order()
{

  addQoI( __SED_QOI_INTEGRATED_FORCE );
  addQoI( __SED_QOI_MIN_DECAY_TIME );

}

//********************************************************************
// modelForcingSimpleExponentialDecay::add_spec_scenarios
//********************************************************************
template <class P_V, class P_M>
  void modelForcingSimpleExponentialDecay<P_V,P_M>::add_spec_scenarios_in_order()
{

  addScn( __SED_SCN_TIME );

}


#endif // __MODEL_FORCING_SIMPLE_EXPONENTIAL_DECAY_H
