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
 * coupled_OscillatorOLS_ForceSLD.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __COUPLED_OSCILLATOR_OLS_FORCE_SLD_H
#define __COUPLED_OSCILLATOR_OLS_FORCE_SLD_H


#include "uqModel.h"
#include "uqDataSet.h"
#include "modelOscillatorLinearSpring.h"
#include "oscillatorLinearSpring.h"
#include "modelForcingSimpleLinearDecay.h"
#include "forcingSimpleLinearDecay.h"


//********************************************************************
// coupled_OscillatorOLS_ForceSLD Definitions
//********************************************************************

// Model name
#define __OLS_SLD_NAME                        "OscOLS_ForceSLD"
// Noise
// #define __OLS_SLD_PARAM_SIGMA                 "param_Osc_OLS_Force_SLD_sigma"
// #define __OLS_SLD_SIGMA                       0.05


//********************************************************************
// coupled_OscillatorOLS_ForceSLD Class
//********************************************************************
template <class P_V, class P_M>
  class coupled_OscillatorOLS_ForceSLD : public uqModel<P_V,P_M>
{

 public:

  coupled_OscillatorOLS_ForceSLD( const uqBaseEnvironmentClass& env );
  coupled_OscillatorOLS_ForceSLD( const uqBaseEnvironmentClass& env, 
			       std::string modelName );
  virtual ~coupled_OscillatorOLS_ForceSLD();

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
// coupled_OscillatorOLS_ForceSLD::coupled_OscillatorOLS_ForceSLD
//********************************************************************
template <class P_V, class P_M>
  coupled_OscillatorOLS_ForceSLD<P_V,P_M>::
  coupled_OscillatorOLS_ForceSLD( const uqBaseEnvironmentClass& env )
  : uqModel<P_V,P_M>( env, __OLS_SLD_NAME )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::coupled_OscillatorOLS_ForceSLD
//********************************************************************
template <class P_V, class P_M>
  coupled_OscillatorOLS_ForceSLD<P_V,P_M>::
  coupled_OscillatorOLS_ForceSLD( const uqBaseEnvironmentClass& env, 
			       std::string modelName )
  : uqModel<P_V,P_M>( env, modelName )
  {
    
    add_spec_params_in_order();
    add_spec_observables_in_order();
    add_spec_qois_in_order();
    add_spec_scenarios_in_order();

  }

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::~coupled_OscillatorOLS_ForceSLD
//********************************************************************
template <class P_V, class P_M>
  coupled_OscillatorOLS_ForceSLD<P_V,P_M>::~coupled_OscillatorOLS_ForceSLD()
{}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::my_forward_problem
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::
  my_forward_problem( P_V& predVec, 
		      std::map<std::string, double> &map_act_params, 
		      std::map<std::string, double> &map_data_point ) const
{

  // add the coupling indicator
  std::map<std::string, double> new_map_act_params( map_act_params );
  new_map_act_params[ __SLD_NAME ] = 1.0;

  // get predVec from Oscillator
  oscillatorLinearSpring_obs<P_V,P_M>( predVec, 
				       map_act_params, 
				       map_data_point );

  // NOTE: one can also augument the observations from both systems
  // in this case only observables from the first systems are used

}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::my_noise_covariance_matrix
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::
  my_noise_covariance_matrix( P_M& covMat, 
			      const P_V& predVec,
			      std::map<std::string, double>& map_act_params ) const
{

  // 1. get it from a file

  // 2. compute it from the parameters
  //  double param_sigma = map_act_params[ __OLS_SLD_PARAM_SIGMA ];
  //  covMat(0,0) = pow( param_sigma, 2.0 );

  // 3. compute it as a % from predicted values

  // 4. use a fixed value
  //  covMat(0,0) = pow( __OLS_SLD_SIGMA, 2.0 );

  UQ_FATAL_TEST_MACRO( true,
		       0,
		       "coupled_OscillatorOLS_ForceSLD::my_noise_covariance_matrix",
		       "implement me" );

}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::my_qoi_model
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::
  my_qoi_model( P_V& qoiPredVec,
		std::map<std::string, double>& map_scn,
		std::map<std::string, double>& map_act_params ) const
{

  // add the coupling indicator
  std::map<std::string, double> new_map_act_params( map_act_params );
  new_map_act_params[ __SLD_NAME ] = 1.0;

  // get the qoi from Oscillator
  uqVectorSpaceClass<P_V,P_M> qoiSpace_oscillator( m_env, 
						   "qoi_oscillator",
						   __OLS_QOI_NUMBER,
						   NULL );
  P_V qoiPredVec_oscillator( qoiSpace_oscillator.zeroVector() );

  oscillatorLinearSpring_qoi<P_V,P_M>( qoiPredVec_oscillator, 
				       map_scn,
				       new_map_act_params );

  // get the qoi from Forcing
  uqVectorSpaceClass<P_V,P_M> qoiSpace_forcing( m_env, 
						"qoi_forcing",
						__SLD_QOI_NUMBER,
						NULL );
  P_V qoiPredVec_forcing( qoiSpace_forcing.zeroVector() );

  simpleLinearDecayForcing_qoi<P_V,P_M>( qoiPredVec_forcing, 
					 map_scn,
					 map_act_params ); 

  // augument the two QoIs vectors
  qoiPredVec.cwSetConcatenated( qoiPredVec_oscillator,
				qoiPredVec_forcing );

}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::add_spec_params
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::add_spec_params_in_order()
{

  // from oscillator OLS
  addParam( __OLS_PARAM_C );
  addParam( __OLS_PARAM_K10 );
  addParam( __OLS_PARAM_SIGMA );

  // from forcing SLD
  addParam( __SLD_PARAM_F0 );
  addParam( __SLD_PARAM_TAU );
  addParam( __SLD_PARAM_SIGMA );

  // specific params
  //  addParam( __OLS_SLD_PARAM_SIGMA );

}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::add_spec_observables_in_order
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::add_spec_observables_in_order()
{

  addObservable( __OLS_OBS_KINETIC_ENERGY, __OLS_OBS_KINETIC_ENERGY_LN );

}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::add_spec_qois_in_order
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::add_spec_qois_in_order()
{

  // from oscillator OLS
  addQoI( __OLS_QOI_TIME_TO_REST );
  addQoI( __OLS_QOI_MAX_DISPLACEMENT );
  addQoI( __OLS_QOI_MAX_VELOCITY );

  // from forcing SLD
  addQoI( __SLD_QOI_INTEGRATED_FORCE );
  addQoI( __SLD_QOI_MIN_DECAY_TIME );

}

//********************************************************************
// coupled_OscillatorOLS_ForceSLD::add_spec_scenarios
//********************************************************************
template <class P_V, class P_M>
  void coupled_OscillatorOLS_ForceSLD<P_V,P_M>::add_spec_scenarios_in_order()
{

  addScn( __OLS_SCN_TIME );

}


#endif // __COUPLED_OSCILLATOR_OLS_FORCE_SLD_H

