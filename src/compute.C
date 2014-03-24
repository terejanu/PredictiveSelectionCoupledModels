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
 * compute.C
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "compute.h"
#include "uqModel.h"
#include "uqParameter.h"
#include "uqDataSet.h"
#include "uqStatsAlg.h"
#include "uqGslVector.h"
#include "uqGslMatrix.h"
#include "uqBagModels.h"
#include <utility>

#include "modelForcingSimpleLinearDecay.h"
#include "modelForcingSimpleExponentialDecay.h"
#include "modelForcingOscillatoryLinearDecay.h"
#include "modelForcingOscillatoryExponentialDecay.h"

#include "modelOscillatorLinearSpring.h"
#include "modelOscillatorCubicSpring.h"
#include "modelOscillatorQuinticSpring.h"

#include "coupled_OscillatorOLS_ForceSLD.h"
#include "coupled_OscillatorOLS_ForceSED.h"
#include "coupled_OscillatorOLS_ForceOLD.h"
#include "coupled_OscillatorOLS_ForceOED.h"

#include "coupled_OscillatorOCS_ForceSLD.h"
#include "coupled_OscillatorOCS_ForceSED.h"
#include "coupled_OscillatorOCS_ForceOLD.h"
#include "coupled_OscillatorOCS_ForceOED.h"

#include "coupled_OscillatorOQS_ForceSLD.h"
#include "coupled_OscillatorOQS_ForceSED.h"
#include "coupled_OscillatorOQS_ForceOLD.h"
#include "coupled_OscillatorOQS_ForceOED.h"

//********************************************************************
// Compute
//********************************************************************
void compute( const uqFullEnvironmentClass& env ) 
{




  //********************************************************************
  //********************************************************************
  //********************************************************************
  // (I) FORCING FUNCTIONS
  //********************************************************************
  //********************************************************************
  //********************************************************************


  //********************************************************************
  // Forcing Data: load the experimental data set
  //********************************************************************
  uqDataSet myDataSet_forcing( "./input_files/forcing_meas.txt" );

  //********************************************************************
  // Forcing Model: Simple Linear Decay (SLD)
  //********************************************************************
  modelForcingSimpleLinearDecay<uqGslVectorClass,uqGslMatrixClass> 
    myForcing_SLD( env );

  myForcing_SLD.setParam( "param_Force_SLD_F0", 0.1, 1.0, false );
  myForcing_SLD.setParam( "param_Force_SLD_tau", PI*5.0, 7.0*PI, false );
  myForcing_SLD.setParam( "param_Force_SLD_sigma", -2.0, 2.0, true );
  myForcing_SLD.finishSetParams();
  myForcing_SLD.printInfo();
  myForcing_SLD.saveInfo();

  //********************************************************************
  // Forcing Model: Simple Exponential Decay (SED)
  //********************************************************************
  modelForcingSimpleExponentialDecay<uqGslVectorClass,uqGslMatrixClass> 
    myForcing_SED( env );

  //  myForcing_SED.setParam( "param_Force_SED_F0", 0.5, 2.0, false );
  //  myForcing_SED.setParam( "param_Force_SED_tau", PI/3.0, 4.0*PI, false );

  myForcing_SED.setParam( "param_Force_SED_F0", 0.1, 2.0, false );
  myForcing_SED.setParam( "param_Force_SED_tau", PI, 10.0*PI, false );

  myForcing_SED.setParam( "param_Force_SED_sigma", -2.0, 2.0, true );
  myForcing_SED.finishSetParams();
  myForcing_SED.printInfo();
  myForcing_SED.saveInfo();

  //********************************************************************
  // Forcing Model: Oscillatory Linear Decay (OLD)
  //********************************************************************
  modelForcingOscillatoryLinearDecay<uqGslVectorClass,uqGslMatrixClass> 
    myForcing_OLD( env );

  myForcing_OLD.setParam( "param_Force_OLD_F0", 0.1, 2.0, false );
  myForcing_OLD.setParam( "param_Force_OLD_tau", PI*5.0, 10.0*PI, false );
  myForcing_OLD.setParam( "param_Force_OLD_alpha", 0.1, 0.4, false );
  myForcing_OLD.setParam( "param_Force_OLD_omega", 1.9, 2.1, false );
  myForcing_OLD.setParam( "param_Force_OLD_sigma", -2.0, 2.0, true );
  myForcing_OLD.finishSetParams();
  myForcing_OLD.printInfo();
  myForcing_OLD.saveInfo();

  //********************************************************************
  // Forcing Model: Oscillatory Exponential Decay (OED)
  //********************************************************************
  modelForcingOscillatoryExponentialDecay<uqGslVectorClass,uqGslMatrixClass> 
    myForcing_OED( env );

  myForcing_OED.setParam( "param_Force_OED_F0", 0.1, 2.0, false );
  myForcing_OED.setParam( "param_Force_OED_tau", PI, 10.0*PI, false );
  myForcing_OED.setParam( "param_Force_OED_alpha", 0.1, 0.4, false );
  myForcing_OED.setParam( "param_Force_OED_omega", 1.9, 2.1, false );
  myForcing_OED.setParam( "param_Force_OED_sigma", -2.0, 2.0, true );
  myForcing_OED.finishSetParams();
  myForcing_OED.printInfo();
  myForcing_OED.saveInfo();

  //********************************************************************
  // Forcing QoI: Prediction scenario & Prediction QoIs 
  //********************************************************************
  qoiDataStruct qoiForcing_DS;

  // NOTE: for integrated quantities <scn_time> is irelevant
  // because we use much larger values
  qoiForcing_DS.m_qoi_scn.push_back( std::make_pair("scn_time", 25.0) ); 

  qoiForcing_DS.m_qois.push_back( "qoi_integrated_force" );
  qoiForcing_DS.m_qois.push_back( "qoi_min_decay_time" );

  //********************************************************************
  // Forcing: Bag of models
  //********************************************************************
  uqBagModels<uqGslVectorClass,uqGslMatrixClass> allModelsForcing( env, "forcing" );

  // add QoI
  allModelsForcing.addQoI( qoiForcing_DS );

  // add Models
  //  allModelsForcing.addModel( &myForcing_SLD );
  allModelsForcing.addModel( &myForcing_SED );
  allModelsForcing.addModel( &myForcing_OLD );
  //  allModelsForcing.addModel( &myForcing_OED );

  //********************************************************************
  // Forcing: Inverse Problems
  //********************************************************************
  for( unsigned int i = 0; i < allModelsForcing.m_bagModels.size(); ++i )
    {
      solveIP( allModelsForcing.m_bagModels[i], myDataSet_forcing );
    }

  //********************************************************************
  // Forcing: Forward Problems
  //********************************************************************
  for( unsigned int i = 0; i < allModelsForcing.m_bagModels.size(); ++i )
    {
      solveFP( allModelsForcing, i );
    }

  //********************************************************************
  // Forcing: Compute evidence & combined predictive distribution
  //********************************************************************
  allModelsForcing.compute_model_plausibility();    // evidence + model plausibility
  allModelsForcing.compute_combined_qoi();          // QoI predictive distribution
  allModelsForcing.save_combined_qoi();             // save chain combined QoI in a file
  allModelsForcing.compute_qoi_aware_ev();          // QoI aware evidence







  //********************************************************************
  //********************************************************************
  //********************************************************************
  // (II) OSCILLATORS
  //********************************************************************
  //********************************************************************
  //********************************************************************



  //********************************************************************
  // Oscillator Data: load the experimental data set
  //********************************************************************
  uqDataSet myDataSet_oscillator( "./input_files/oscillator_meas.txt" );

  //********************************************************************
  // Oscillator Model: Oscillator Linear Spring (OLS)
  //********************************************************************
  modelOscillatorLinearSpring<uqGslVectorClass,uqGslMatrixClass> 
    myOscillator_OLS( env );

  myOscillator_OLS.setParam( "param_Osc_OLS_c", 0.05, 0.15, false );
  myOscillator_OLS.setParam( "param_Osc_OLS_k10", 3.0, 5.0, false );
  myOscillator_OLS.setParam( "param_Osc_OLS_sigma", -3.0, 1.0, true );
  myOscillator_OLS.finishSetParams();
  myOscillator_OLS.printInfo();
  myOscillator_OLS.saveInfo();

  //********************************************************************
  // Oscillator Model: Oscillator Cubic Spring (OCS)
  //********************************************************************
  modelOscillatorCubicSpring<uqGslVectorClass,uqGslMatrixClass> 
    myOscillator_OCS( env );

  myOscillator_OCS.setParam( "param_Osc_OCS_c", 0.05, 0.15, false );
  myOscillator_OCS.setParam( "param_Osc_OCS_k30", 3.0, 5.0, false );
  myOscillator_OCS.setParam( "param_Osc_OCS_k32", -1.0, 1.0, false );
  myOscillator_OCS.setParam( "param_Osc_OCS_sigma", -3.0, 1.0, true );
  myOscillator_OCS.finishSetParams();
  myOscillator_OCS.printInfo();
  myOscillator_OCS.saveInfo();

  //********************************************************************
  // Oscillator Model: Oscillator Quintic Spring (OQS)
  //********************************************************************
  modelOscillatorQuinticSpring<uqGslVectorClass,uqGslMatrixClass> 
    myOscillator_OQS( env );

  myOscillator_OQS.setParam( "param_Osc_OQS_c", 0.05, 0.15, false );
  myOscillator_OQS.setParam( "param_Osc_OQS_k50", 3.0, 5.0, false );
  myOscillator_OQS.setParam( "param_Osc_OQS_k52", -6.0, -4.0, false );
  myOscillator_OQS.setParam( "param_Osc_OQS_k54", 0.5, 1.5, false );
  myOscillator_OQS.setParam( "param_Osc_OQS_sigma", -3.0, 1.0, true );
  myOscillator_OQS.finishSetParams();
  myOscillator_OQS.printInfo();
  myOscillator_OQS.saveInfo();

  //********************************************************************
  // Oscillator QoI: Prediction scenario & Prediction QoIs 
  //********************************************************************
  qoiDataStruct qoiOscillator_DS;

  // NOTE: for integrated quantities <scn_time> is irelevant
  // because we use much larger values
  qoiOscillator_DS.m_qoi_scn.push_back( std::make_pair("scn_time", 25.0) ); 

  qoiOscillator_DS.m_qois.push_back( "qoi_time_to_rest" );
  qoiOscillator_DS.m_qois.push_back( "qoi_max_displacement" );
  qoiOscillator_DS.m_qois.push_back( "qoi_max_velocity" );

  //********************************************************************
  // Oscillator: Bag of models
  //********************************************************************
  uqBagModels<uqGslVectorClass,uqGslMatrixClass> allModelsOscillators( env, "oscillator" );

  // add QoI
  allModelsOscillators.addQoI( qoiOscillator_DS );

  // add Models
  allModelsOscillators.addModel( &myOscillator_OLS );
  allModelsOscillators.addModel( &myOscillator_OCS );
  //  allModelsOscillators.addModel( &myOscillator_OQS );

  //********************************************************************
  // Oscillator: Inverse Problems
  //********************************************************************
  for( unsigned int i = 0; i < allModelsOscillators.m_bagModels.size(); ++i )
    {
      solveIP( allModelsOscillators.m_bagModels[i], myDataSet_oscillator );
    }

  //********************************************************************
  // Oscillator: Forward Problems
  //********************************************************************
  for( unsigned int i = 0; i < allModelsOscillators.m_bagModels.size(); ++i )
    {
      solveFP( allModelsOscillators, i );
    }

  //********************************************************************
  // Oscillator: Compute evidence & combined predictive distribution
  //********************************************************************
  allModelsOscillators.compute_model_plausibility();  // evidence + model plausibility
  allModelsOscillators.compute_combined_qoi();        // QoI predictive distribution
  allModelsOscillators.save_combined_qoi();           // save chain combined QoI in a file
  allModelsOscillators.compute_qoi_aware_ev();        // QoI aware evidence





  //********************************************************************
  //********************************************************************
  //********************************************************************
  // (III) COUPLED MODELS: OSCILLATORS + FORCING FUNCTIONS
  //********************************************************************
  //********************************************************************
  //********************************************************************

  // --------- OLS

  //********************************************************************
  // Coupled Model: Oscillator (OLS) - Forcing (SLD)
  //********************************************************************
  coupled_OscillatorOLS_ForceSLD<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OLS_SLD( env );

  // set parameters
  myCoupled_OLS_SLD.setParamsFromModel( myForcing_SLD );
  myCoupled_OLS_SLD.setParamsFromModel( myOscillator_OLS );
  //  myCoupled_OLS_SLD.setParam( "param_Osc_OLS_Force_SLD_sigma", -3.0, 1.0, true );
  myCoupled_OLS_SLD.finishSetParams();
  myCoupled_OLS_SLD.printInfo();
  myCoupled_OLS_SLD.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OLS) - Forcing (SED)
  //********************************************************************
  coupled_OscillatorOLS_ForceSED<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OLS_SED( env );

  // set parameters
  myCoupled_OLS_SED.setParamsFromModel( myForcing_SED );
  myCoupled_OLS_SED.setParamsFromModel( myOscillator_OLS );
  //  myCoupled_OLS_SED.setParam( "param_Osc_OLS_Force_SED_sigma", -3.0, 1.0, true );
  myCoupled_OLS_SED.finishSetParams();
  myCoupled_OLS_SED.printInfo();
  myCoupled_OLS_SED.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OLS) - Forcing (OLD)
  //********************************************************************
  coupled_OscillatorOLS_ForceOLD<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OLS_OLD( env );

  // set parameters
  myCoupled_OLS_OLD.setParamsFromModel( myForcing_OLD );
  myCoupled_OLS_OLD.setParamsFromModel( myOscillator_OLS );
  //  myCoupled_OLS_OLD.setParam( "param_Osc_OLS_Force_OLD_sigma", -3.0, 1.0, true );
  myCoupled_OLS_OLD.finishSetParams();
  myCoupled_OLS_OLD.printInfo();
  myCoupled_OLS_OLD.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OLS) - Forcing (OED)
  //********************************************************************
  coupled_OscillatorOLS_ForceOED<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OLS_OED( env );

  // set parameters
  myCoupled_OLS_OED.setParamsFromModel( myForcing_OED );
  myCoupled_OLS_OED.setParamsFromModel( myOscillator_OLS );
  //  myCoupled_OLS_OED.setParam( "param_Osc_OLS_Force_OED_sigma", -3.0, 1.0, true );
  myCoupled_OLS_OED.finishSetParams();
  myCoupled_OLS_OED.printInfo();
  myCoupled_OLS_OED.saveInfo();

  // --------- OCS

  //********************************************************************
  // Coupled Model: Oscillator (OCS) - Forcing (SLD)
  //********************************************************************
  coupled_OscillatorOCS_ForceSLD<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OCS_SLD( env );

  // set parameters
  myCoupled_OCS_SLD.setParamsFromModel( myForcing_SLD );
  myCoupled_OCS_SLD.setParamsFromModel( myOscillator_OCS );
  //  myCoupled_OCS_SLD.setParam( "param_Osc_OCS_Force_SLD_sigma", -3.0, 1.0, true );
  myCoupled_OCS_SLD.finishSetParams();
  myCoupled_OCS_SLD.printInfo();
  myCoupled_OCS_SLD.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OCS) - Forcing (SED)
  //********************************************************************
  coupled_OscillatorOCS_ForceSED<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OCS_SED( env );

  // set parameters
  myCoupled_OCS_SED.setParamsFromModel( myForcing_SED );
  myCoupled_OCS_SED.setParamsFromModel( myOscillator_OCS );
  //  myCoupled_OCS_SED.setParam( "param_Osc_OCS_Force_SED_sigma", -3.0, 1.0, true );
  myCoupled_OCS_SED.finishSetParams();
  myCoupled_OCS_SED.printInfo();
  myCoupled_OCS_SED.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OCS) - Forcing (OLD)
  //********************************************************************
  coupled_OscillatorOCS_ForceOLD<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OCS_OLD( env );

  // set parameters
  myCoupled_OCS_OLD.setParamsFromModel( myForcing_OLD );
  myCoupled_OCS_OLD.setParamsFromModel( myOscillator_OCS );
  //  myCoupled_OCS_OLD.setParam( "param_Osc_OCS_Force_OLD_sigma", -3.0, 1.0, true );
  myCoupled_OCS_OLD.finishSetParams();
  myCoupled_OCS_OLD.printInfo();
  myCoupled_OCS_OLD.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OCS) - Forcing (OED)
  //********************************************************************
  coupled_OscillatorOCS_ForceOED<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OCS_OED( env );

  // set parameters
  myCoupled_OCS_OED.setParamsFromModel( myForcing_OED );
  myCoupled_OCS_OED.setParamsFromModel( myOscillator_OCS );
  //  myCoupled_OCS_OED.setParam( "param_Osc_OCS_Force_OED_sigma", -3.0, 1.0, true );
  myCoupled_OCS_OED.finishSetParams();
  myCoupled_OCS_OED.printInfo();
  myCoupled_OCS_OED.saveInfo();


  // --------- OQS

  //********************************************************************
  // Coupled Model: Oscillator (OQS) - Forcing (SLD)
  //********************************************************************
  coupled_OscillatorOQS_ForceSLD<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OQS_SLD( env );

  // set parameters
  myCoupled_OQS_SLD.setParamsFromModel( myForcing_SLD );
  myCoupled_OQS_SLD.setParamsFromModel( myOscillator_OQS );
  //  myCoupled_OQS_SLD.setParam( "param_Osc_OQS_Force_SLD_sigma", -3.0, 1.0, true );
  myCoupled_OQS_SLD.finishSetParams();
  myCoupled_OQS_SLD.printInfo();
  myCoupled_OQS_SLD.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OQS) - Forcing (SED)
  //********************************************************************
  coupled_OscillatorOQS_ForceSED<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OQS_SED( env );

  // set parameters
  myCoupled_OQS_SED.setParamsFromModel( myForcing_SED );
  myCoupled_OQS_SED.setParamsFromModel( myOscillator_OQS );
  //  myCoupled_OQS_SED.setParam( "param_Osc_OQS_Force_SED_sigma", -3.0, 1.0, true );
  myCoupled_OQS_SED.finishSetParams();
  myCoupled_OQS_SED.printInfo();
  myCoupled_OQS_SED.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OQS) - Forcing (OLD)
  //********************************************************************
  coupled_OscillatorOQS_ForceOLD<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OQS_OLD( env );

  // set parameters
  myCoupled_OQS_OLD.setParamsFromModel( myForcing_OLD );
  myCoupled_OQS_OLD.setParamsFromModel( myOscillator_OQS );
  //  myCoupled_OQS_OLD.setParam( "param_Osc_OQS_Force_OLD_sigma", -3.0, 1.0, true );
  myCoupled_OQS_OLD.finishSetParams();
  myCoupled_OQS_OLD.printInfo();
  myCoupled_OQS_OLD.saveInfo();

  //********************************************************************
  // Coupled Model: Oscillator (OQS) - Forcing (OED)
  //********************************************************************
  coupled_OscillatorOQS_ForceOED<uqGslVectorClass,uqGslMatrixClass> 
    myCoupled_OQS_OED( env );

  // set parameters
  myCoupled_OQS_OED.setParamsFromModel( myForcing_OED );
  myCoupled_OQS_OED.setParamsFromModel( myOscillator_OQS );
  //  myCoupled_OQS_OED.setParam( "param_Osc_OQS_Force_OED_sigma", -3.0, 1.0, true );
  myCoupled_OQS_OED.finishSetParams();
  myCoupled_OQS_OED.printInfo();
  myCoupled_OQS_OED.saveInfo();


  //********************************************************************
  // Coupled Model: Oscillator QoI - Prediction scenario & Prediction QoIs 
  //********************************************************************
  qoiDataStruct qoiCoupled_DS;

  // NOTE: for integrated quantities <scn_time> is irelevant
  // because we use much larger values
  qoiCoupled_DS.m_qoi_scn.push_back( std::make_pair("scn_time", 25.0) ); 

  qoiCoupled_DS.m_qois.push_back( "qoi_time_to_rest" );
  qoiCoupled_DS.m_qois.push_back( "qoi_max_displacement" );
  qoiCoupled_DS.m_qois.push_back( "qoi_max_velocity" );
  qoiCoupled_DS.m_qois.push_back( "qoi_integrated_force" );
  qoiCoupled_DS.m_qois.push_back( "qoi_min_decay_time" );

  //********************************************************************
  // Coupled Model: Bag of models
  //********************************************************************
  uqBagModels<uqGslVectorClass,uqGslMatrixClass> allModelsCoupled( env, "coupled" );

  // add QoI
  allModelsCoupled.addQoI( qoiCoupled_DS );

  // add Models
  //  allModelsCoupled.addModel( &myCoupled_OLS_SLD );
  allModelsCoupled.addModel( &myCoupled_OLS_SED );
  allModelsCoupled.addModel( &myCoupled_OLS_OLD );
  //  allModelsCoupled.addModel( &myCoupled_OLS_OED );

  //  allModelsCoupled.addModel( &myCoupled_OCS_SLD );
  allModelsCoupled.addModel( &myCoupled_OCS_SED );
  allModelsCoupled.addModel( &myCoupled_OCS_OLD );
  //  allModelsCoupled.addModel( &myCoupled_OCS_OED );

  //  allModelsCoupled.addModel( &myCoupled_OQS_SLD );
  //  allModelsCoupled.addModel( &myCoupled_OQS_SED );
  //  allModelsCoupled.addModel( &myCoupled_OQS_OLD );
  //  allModelsCoupled.addModel( &myCoupled_OQS_OED );



  //********************************************************************
  // Coupled Model: Construct PostRV & Model Plausibilities
  //********************************************************************

  // ------- OLS

  /*
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OLS_SLD.getName(),
			  allModelsOscillators,  myOscillator_OLS.getName(),
			  allModelsForcing,      myForcing_SLD.getName() );
  */
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OLS_SED.getName(),
			  allModelsOscillators,  myOscillator_OLS.getName(),
			  allModelsForcing,      myForcing_SED.getName() );

  constructCoupledPostRV( allModelsCoupled,      myCoupled_OLS_OLD.getName(),
			  allModelsOscillators,  myOscillator_OLS.getName(),
			  allModelsForcing,      myForcing_OLD.getName() );
  /*
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OLS_OED.getName(),
			  allModelsOscillators,  myOscillator_OLS.getName(),
			  allModelsForcing,      myForcing_OED.getName() );
  */
  // ------- OCS
  /*
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OCS_SLD.getName(),
			  allModelsOscillators,  myOscillator_OCS.getName(),
			  allModelsForcing,      myForcing_SLD.getName() );
  */
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OCS_SED.getName(),
			  allModelsOscillators,  myOscillator_OCS.getName(),
			  allModelsForcing,      myForcing_SED.getName() );

  constructCoupledPostRV( allModelsCoupled,      myCoupled_OCS_OLD.getName(),
			  allModelsOscillators,  myOscillator_OCS.getName(),
			  allModelsForcing,      myForcing_OLD.getName() );
  /*
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OCS_OED.getName(),
			  allModelsOscillators,  myOscillator_OCS.getName(),
			  allModelsForcing,      myForcing_OED.getName() );
  */
  // ------- OQS
  /*
  constructCoupledPostRV( allModelsCoupled,      myCoupled_OQS_SLD.getName(),
			  allModelsOscillators,  myOscillator_OQS.getName(),
			  allModelsForcing,      myForcing_SLD.getName() );

  constructCoupledPostRV( allModelsCoupled,      myCoupled_OQS_SED.getName(),
			  allModelsOscillators,  myOscillator_OQS.getName(),
			  allModelsForcing,      myForcing_SED.getName() );

  constructCoupledPostRV( allModelsCoupled,      myCoupled_OQS_OLD.getName(),
			  allModelsOscillators,  myOscillator_OQS.getName(),
			  allModelsForcing,      myForcing_OLD.getName() );

  constructCoupledPostRV( allModelsCoupled,      myCoupled_OQS_OED.getName(),
			  allModelsOscillators,  myOscillator_OQS.getName(),
			  allModelsForcing,      myForcing_OED.getName() );
  */
  //********************************************************************
  // Coupled Model: Forward Problems
  //********************************************************************
  for( unsigned int i = 0; i < allModelsCoupled.m_bagModels.size(); ++i )
    {
      solveFP( allModelsCoupled, i );
    }

  //********************************************************************
  // Coupled Model: Compute Conditional Mutual Information
  //********************************************************************

  std::vector<std::string> surrogate_qois;
  surrogate_qois.push_back( "qoi_integrated_force" );
  surrogate_qois.push_back( "qoi_min_decay_time" );

  allModelsCoupled.compute_conditional_mi( "qoi_time_to_rest",
					   surrogate_qois );


  //********************************************************************
  // Coupled Model: Compute evidence & combined predictive distribution
  //********************************************************************

  allModelsCoupled.compute_combined_qoi();          // QoI predictive distribution
  allModelsCoupled.save_combined_qoi();             // save chain combined QoI in a file
  allModelsCoupled.compute_qoi_aware_ev();          // QoI aware evidence






  //********************************************************************
  //********************************************************************
  //********************************************************************
  // (IV) SUMMARY RESULTS
  //********************************************************************
  //********************************************************************
  //********************************************************************


  //********************************************************************
  // Print evidence
  //********************************************************************
  allModelsForcing.print_evidences();
  allModelsOscillators.print_evidences();
  allModelsCoupled.print_evidences();
  allModelsCoupled.print_conditional_mi();


}


