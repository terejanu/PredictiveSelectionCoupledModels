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
 * uqModel.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MODEL_H
#define __UQ_MODEL_H

#include <uqStatisticalInverseProblem.h>
#include "uqParameter.h"
#include "uqVectorRV.h"
#include "uqGslMatrix.h"
#include <map>
#include "uqDataSet.h"

#ifndef PI
#define PI 3.1415926535897932
#endif // PI

//********************************************************************
// Data Structures
//********************************************************************
struct obsStruct
{
  std::string   obsName;
  bool          obsLn;
};

struct qoiDataStruct
{
 public:
  std::vector< std::string >                       m_qois;
  std::vector< std::pair<std::string,double> >     m_qoi_scn;
};

//********************************************************************
// uqModel Class
//********************************************************************
template <class P_V, class P_M>
class uqModel
{

 public:

  P_V*                                           paramMins;
  P_V*                                           paramMaxs;
  uqVectorSpaceClass<P_V,P_M>*                   paramSpace;
  uqBoxSubsetClass<P_V,P_M>*                     paramDomain;
  uqUniformVectorRVClass<P_V,P_M>*               priorRV;

  // model
  uqModel( const uqBaseEnvironmentClass& env );
  uqModel( const uqBaseEnvironmentClass& env, std::string modelName );
  virtual ~uqModel();

  // parameters 
  void setParam( std::string name, double min_val, double max_val, bool log_scale );
  bool getParam( std::string name, double &min_val, double &max_val, bool &log_scale ) const;
  void setParamsFromModel( const uqModel& model );
  void finishSetParams();

  virtual void printInfo();
  virtual void saveInfo();
  unsigned int getNumberParams() const;
  unsigned int getNumberObservables() const;
  std::string getName() const;

  void getActualParameters( const P_V& paramValues, std::map<std::string,double>& map_act_params ) const;

  // statistics
  virtual double likelihoodFunction(
				    const P_V&   paramValues,
				    const P_V*   paramDirection,
				    const void*  functionDataPtr,
				    P_V*         gradVector,
				    P_M*         hessianMatrix,
				    P_V*         hessianEffect ) const;
  
  virtual void qoiFunction(
			   const P_V&                paramValues,
			   const P_V*                paramDirection,
			   const void*               functionDataPtr,
			   P_V&                      qoiValues,
			   uqDistArrayClass<P_V*>*   gradVector,
			   uqDistArrayClass<P_M*>*   hessianMatrix,
			   uqDistArrayClass<P_V*>*   hessianEffect) const;

  // test functions
  void test_printActualValue( const P_V& paramValues );  

 protected:

  std::string                                    m_prefix;
  const uqBaseEnvironmentClass&                  m_env;

  unsigned int                                   m_num_params;
  unsigned int                                   m_num_params_set;
  bool                                           m_set_params;
  std::map<std::string,uqParameter>              m_all_params;
  std::vector<std::string>                       m_order_params;
  unsigned int                                   m_num_observables;
  std::vector<obsStruct>                         m_order_observables;
  unsigned int                                   m_num_qois;
  std::vector<std::string>                       m_order_qois;
  unsigned int                                   m_num_scns;
  std::vector<std::string>                       m_order_scns;

  virtual void add_spec_params_in_order() = 0;

  void addParam( std::string name );
  void addObservable( std::string name, bool obsLn );
  void addQoI( std::string name );
  void addScn( std::string name );

  void compute_obs_pred_diff( P_V& diffVec, const P_V& predVec,
			      std::map<std::string, double>& map_data_point ) const;
  double get_correction_term( std::map<std::string, double>& map_data_point ) const;

  void get_act_val_param( const P_V& paramValues, P_V& actualValues );
  virtual double evaluate_loglikelihood( const uqDataSet& myDataSet, 
					 std::map<std::string, double>& map_act_params) const;

  virtual void my_forward_problem( P_V& predVec, 
				   std::map<std::string, double>& map_act_params, 
				   std::map<std::string, double>& map_data_point ) const = 0;
  virtual void my_noise_covariance_matrix( P_M& covMat, 
					   const P_V& predVec,
					   std::map<std::string, double>& map_act_params ) const = 0;

  virtual void my_qoi_model( P_V& qoiPredVec,
			     std::map<std::string, double>& map_scn,
			     std::map<std::string, double>& map_act_params ) const = 0;

};

//********************************************************************
// uqModel::uqModel
//********************************************************************
template <class P_V, class P_M>
  uqModel<P_V,P_M>::uqModel( const uqBaseEnvironmentClass& env )
  : m_prefix( "NONE" ),
    m_env( env ),
    m_num_params( 0 ),
    m_num_params_set( 0 ),
    m_set_params( false ),
    m_num_observables( 0 ),
    m_num_qois( 0 ),
    m_num_scns( 0 )
  {}

//********************************************************************
// uqModel::uqModel
//********************************************************************
template <class P_V, class P_M>
  uqModel<P_V,P_M>::uqModel( const uqBaseEnvironmentClass& env, std::string modelName )
  : m_prefix( modelName ),
    m_env( env ),
    m_num_params( 0 ),
    m_num_params_set( 0 ),
    m_set_params( false ),
    m_num_observables( 0 ),
    m_num_qois( 0 ),
    m_num_scns( 0 )
  {}

//********************************************************************
// uqModel::~uqModel
//********************************************************************
template <class P_V, class P_M>
  uqModel<P_V,P_M>::~uqModel()
{
  if( m_set_params )
    {
      //      delete postRV;
      delete priorRV;
      delete paramDomain;
      delete paramMaxs;
      delete paramMins;
      delete paramSpace;
    }
}

//********************************************************************
// uqModel::likelihoodFunction
//********************************************************************
template <class P_V, class P_M>
  double uqModel<P_V,P_M>::likelihoodFunction(
					       const P_V&   paramValues,
					       const P_V*   paramDirection,
					       const void*  functionDataPtr,
					       P_V*         gradVector,
					       P_M*         hessianMatrix,
					       P_V*         hessianEffect ) const
{
  // my data for the likelihood function
  std::map<std::string, double> map_act_params;

  // get the data set
  const uqDataSet& myDataSet = *((uqDataSet *) functionDataPtr);

  // get the parameters
  getActualParameters( paramValues, map_act_params );

  // compute the likelihood
  double loglike_val = this->evaluate_loglikelihood( myDataSet, map_act_params );

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  return loglike_val;

#else
  return -2 * loglike_val;

#endif
}

//********************************************************************
// uqModel::qoiFunction
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::qoiFunction( const P_V&                paramValues,
				      const P_V*                paramDirection,
				      const void*               functionDataPtr,
				      P_V&                      qoiValues,
				      uqDistArrayClass<P_V*>*   gradVector,
				      uqDistArrayClass<P_M*>*   hessianMatrix,
				      uqDistArrayClass<P_V*>*   hessianEffect) const
{

  // my data for the likelihood function
  std::map<std::string, double> map_act_params;

  // get the parameters
  getActualParameters( paramValues, map_act_params );

  // my scenario
  qoiDataStruct qDS = *((qoiDataStruct *) functionDataPtr);
  const std::vector< std::string >&  m_qois = qDS.m_qois;
  const std::vector< std::pair<std::string,double> >& m_qoi_scn = qDS.m_qoi_scn;

  // helper variables
  std::vector< unsigned int > corresp_vec;
  std::map<std::string, double> map_scn;
  
  // check that all the qoi's requested can be returned
  for( unsigned int i = 0; i < m_qois.size(); ++i )
    {
      bool qoi_found = false;

      for( unsigned int j = 0; (j < m_order_qois.size()) && !qoi_found ; ++j )
	{
	  if( m_qois[i] == m_order_qois[j] )
	    {
	      qoi_found = true;
	      corresp_vec.push_back( j );
	    }
	}

      UQ_FATAL_TEST_MACRO( !qoi_found,
			   0,
			   "uqModel::qoiFunction",
			   "One of the qois requested is not present in the model" );
    }

  // create the qoi map
  for( unsigned int i = 0; i < m_qoi_scn.size(); ++i )
    {
      map_scn[ m_qoi_scn[i].first ] = m_qoi_scn[i].second;
    }

  // compute the model qois
  uqVectorSpaceClass<P_V,P_M>
    qoiSpace(m_env, "qoi_", m_num_qois, NULL);
  P_V qoiPredVec( qoiSpace.zeroVector() );

  this->my_qoi_model( qoiPredVec, map_scn, map_act_params );

  // arrage the response value
  for( unsigned int i = 0; i < m_qois.size(); ++i )
    {
      qoiValues[ i ] = qoiPredVec[ corresp_vec[i] ];
    }  

}

//********************************************************************
// uqLinForcing::evaluate_likelihood
//********************************************************************
template <class P_V, class P_M>
  double uqModel<P_V,P_M>::evaluate_loglikelihood( const uqDataSet& myDataSet, 
						   std::map<std::string, double>& map_act_params ) const
{

  double loglike_val = 0;
  std::map<std::string, double> map_data_point;
  unsigned int no_observables = getNumberObservables();

  uqVectorSpaceClass<P_V,P_M> obsSpace( m_env, 
					("obs_" + m_prefix).c_str(),
					no_observables,
					NULL );

  P_V predVec( obsSpace.zeroVector() );
  P_V diffVec( obsSpace.zeroVector() );
  P_M* covMat = obsSpace.newMatrix();

  for( unsigned int i = 0; i < myDataSet.getDataSetSize(); ++i )
    {
      // get a data point
      myDataSet.getMapDataPoint( i, map_data_point );

      // compute forward problem
      this->my_forward_problem( predVec, map_act_params, map_data_point );

      // get the covariance matrix
      this->my_noise_covariance_matrix( *covMat, predVec, map_act_params );

      // compute the error
      compute_obs_pred_diff( diffVec, predVec, map_data_point );

      // get correction term due to potential multiuplicative noise
      double corr_term = get_correction_term( map_data_point );

      //-----------------------------------------------------------------
      /*
      std::cout << std::endl << "passed " << i << " ..." << std::endl;
      std::cout << "param values : " << std::endl;
      for( std::map<std::string,double>::iterator it = map_act_params.begin();
	   it != map_act_params.end(); ++it )
	{
	  std::cout << "\t" << it->first << " = " << it->second << std::endl;
	}
      std::cout << "data point : " << std::endl;
      for( std::map<std::string,double>::iterator it = map_data_point.begin();
	   it != map_data_point.end(); ++it )
	{
	  std::cout << "\t" << it->first << " = " << it->second << std::endl;
	}
      std::cout << "pred vec = " << predVec << std::endl;
      std::cout << "diff vec = " << diffVec << std::endl;
      std::cout << "scalar prod = " << scalarProduct( diffVec, covMat->invertMultiply( diffVec ) ) << std::endl;
      std::cout << "cov mat = " << (*covMat) << std::endl;
      std::cout << "lnDet = " << covMat->lnDeterminant() << std::endl;
      std::cout << "corr term = " << corr_term << std::endl;
      */
      //-----------------------------------------------------------------

      // compute likelihood

      double quad_term = scalarProduct( diffVec, covMat->invertMultiply( diffVec ) );
      double cov_lnDet = covMat->lnDeterminant();

      loglike_val += 
	- 0.5 * (double)no_observables * std::log( 2.0 * PI ) * cov_lnDet
	- 0.5 * quad_term - corr_term;

      //-----------------------------------------------------------------
      // std::cout << "loglike_val = "  << loglike_val << std::endl;
      // UQ_FATAL_TEST_MACRO( true, 0, "", "sanity check" );
      //-----------------------------------------------------------------

    }

  // release memory
  delete covMat;

  return loglike_val;

}


//********************************************************************
// uqModel::addParam
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::addParam( std::string name )
{

  UQ_FATAL_TEST_MACRO( m_set_params,
		       0,
		       "uqModel::addParam",
		       "m_set_params == true - Cannot add a new param after all have been set" );

  m_num_params++;
  m_all_params[ name ].setName( name );
  m_order_params.push_back( name );  

}

//********************************************************************
// uqModel::setParam
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::setParam( std::string name, 
				   double min_val, 
				   double max_val, 
				   bool log_scale )
{

  UQ_FATAL_TEST_MACRO( m_set_params,
		       0,
		       "uqModel::setParam",
		       "m_set_params == true - Cannot set another param after all have been set" );

  UQ_FATAL_TEST_MACRO( m_all_params.find( name ) == m_all_params.end(),
		       0,
		       "uqModel::setParam",
		       "Trying to set a parameter that is not in the hashmap" );

  m_all_params[ name ].setParam( min_val, max_val, log_scale );
  m_num_params_set++;

}

//********************************************************************
// uqModel::getParam
//********************************************************************
template <class P_V, class P_M>
  bool uqModel<P_V,P_M>::getParam( std::string name, 
				  double &min_val, 
				  double &max_val, 
				  bool &log_scale ) const
{

  std::map<std::string,uqParameter>::const_iterator it;
  it = m_all_params.find( name );

  if( it == m_all_params.end() )
    {
      return false;    // WARNING - I do not have this parameter
    }
  else
    {
      it->second.getParam( min_val, max_val, log_scale );
      return true;    // SUCCESS - The parameter is in my list
    }

}


//********************************************************************
// uqModel::setParamsFromModel
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::setParamsFromModel( const uqModel& model )
{

  double min_val, max_val;
  bool log_scale, ret_success;

  for( unsigned int i = 0; i < m_num_params; ++i )
    {

      ret_success = model.getParam( m_order_params[i],
				    min_val, max_val, log_scale );
      if( ret_success )
	{
	  this->setParam( m_order_params[i],
			  min_val, max_val, log_scale );
	}

    }

}


//********************************************************************
// uqModel::addObservable
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::addObservable( std::string name, bool obsLn )
{

  obsStruct tmp;
  
  tmp.obsName = name;
  tmp.obsLn = obsLn;

  m_order_observables.push_back( tmp );
  m_num_observables++;

}

//********************************************************************
// uqModel::addQoI
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::addQoI( std::string name )
{

  m_order_qois.push_back( name );
  m_num_qois++;

}

//********************************************************************
// uqModel::addScn
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::addScn( std::string name )
{

  m_order_scns.push_back( name );
  m_num_scns++;

}


//********************************************************************
// uqModel::finishSetParams
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::finishSetParams()
{

  if(!m_set_params)
    {

      UQ_FATAL_TEST_MACRO( m_num_params == 0,
			   0,
			   "uqModel::finishSetParams",
			   "m_num_params == 0 - There are no parameters added in the hashmap" );

      UQ_FATAL_TEST_MACRO( m_num_params != m_num_params_set,
			   0,
			   "uqModel::finishSetParams",
			   "m_num_params != m_num_params_set - Not all the parameters have been set" );

      m_set_params = true;

      paramSpace = new uqVectorSpaceClass<P_V,P_M>( m_env, 
						    ("param_" + m_prefix).c_str(),
						    m_num_params,
						    NULL );
      
      paramMins = new P_V( paramSpace->zeroVector() );
      paramMaxs = new P_V( paramSpace->zeroVector() );
      
      for( unsigned int i = 0; i < m_num_params; ++i )
	{
	  (*paramMins)[ i ] = m_all_params[ m_order_params[ i ] ].getNormMin();
	  (*paramMaxs)[ i ] = m_all_params[ m_order_params[ i ] ].getNormMax();
	}
      
      paramDomain = new uqBoxSubsetClass<P_V,P_M>( ("param_" + m_prefix).c_str(),
						   *paramSpace,
						   *paramMins,
						   *paramMaxs );
      
      priorRV = new uqUniformVectorRVClass<P_V,P_M>( ("prior_" + m_prefix).c_str(),
						     *paramDomain );

    }

}

//********************************************************************
// uqModel::printInfo
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::printInfo()
{

  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "Model name          : " << m_prefix << std::endl;
  std::cout << "Number params       : " << m_num_params << std::endl;
  std::cout << "Number params set   : " << m_num_params_set << std::endl;
  std::cout << "All params are set? : " << m_set_params << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;

  std::cout << "All params      : ";

  for( std::map<std::string,uqParameter>::iterator it = m_all_params.begin();
       it != m_all_params.end(); it++ )
    {
      std::cout << (*it).first << " | ";
    }

  std::cout << std::endl;

  std::cout << "All observables : ";

  for( unsigned int i = 0; i < m_num_observables; ++i )
    {
      std::cout << m_order_observables[i].obsName << " | ";
    }

  std::cout << std::endl;

  std::cout << "All qois        : ";

  for( unsigned int i = 0; i < m_num_qois; ++i ) 
    {
      std::cout << m_order_qois[i] << " | ";
    }

  std::cout << std::endl;

  std::cout << "All scenarios   : ";

  for( unsigned int i = 0; i < m_num_scns; ++i ) 
    {
      std::cout << m_order_scns[i] << " | ";
    }

  std::cout << std::endl;

  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "Parameters set: " << std::endl;

  for( unsigned int i = 0; i < m_num_params_set; ++i )
    {
      m_all_params[ m_order_params[ i ] ].printInfo();
    }
  std::cout << "-----------------------------------------------------" << std::endl;

}


//********************************************************************
// uqModel::saveInfo
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::saveInfo()
{

  std::string fileName = "./outputData/" + m_prefix + "_info.m";

  std::ofstream qoiFile;
  qoiFile.open ( fileName.c_str() );

  qoiFile << "params_min_max_log = zeros(" << m_num_params << ", 3 );" << std::endl;
  qoiFile << "params_min_max_log = [";

  for( unsigned int i = 0; i < m_num_params; ++i )
    {

      double min_val, max_val;
      bool log_scale;

      m_all_params[ m_order_params[ i ] ].getParam( min_val, max_val, log_scale );

      qoiFile << min_val << ", " << max_val << ", ";
      if( log_scale )
	qoiFile << "1";
      else
	qoiFile << "0";

      if( i < (m_num_params-1) )
	qoiFile << ";" << std::endl;
      else
	qoiFile << std::endl;

    }
  
  qoiFile << "];" << std::endl;

  qoiFile.close();

}



//********************************************************************
// uqModel::getNumberParams
//********************************************************************
template <class P_V, class P_M>
  unsigned int uqModel<P_V,P_M>::getNumberParams() const
{
  return m_num_params;
}

//********************************************************************
// uqModel::getNumberObservables
//********************************************************************
template <class P_V, class P_M>
  unsigned int uqModel<P_V,P_M>::getNumberObservables() const
{
  return m_num_observables;
}

//********************************************************************
// uqModel::getName
//********************************************************************
template <class P_V, class P_M>
  std::string uqModel<P_V,P_M>::getName( ) const
{

  return m_prefix;

}

//********************************************************************
// uqModel::getActualParameters
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::getActualParameters( const P_V& paramValues, std::map<std::string,double>& map_act_params ) const
{

  map_act_params.clear();

  std::map<std::string, uqParameter>::const_iterator it;

  UQ_FATAL_TEST_MACRO( !m_set_params,
		       0,
		       "uqModel::getActualParameters",
		       "!m_set_params - Computing actual values is not possible if all the parameters are not set" );
  
  for( unsigned int i = 0; i < m_num_params; ++i )
    {

      it = m_all_params.find( m_order_params[ i ] );
      map_act_params[ m_order_params[ i ] ] = it->second.getActualValue( paramValues[ i ] );
    }

}

//********************************************************************
// uqModel::compute_obs_pred_diff
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::compute_obs_pred_diff( P_V& diffVec, const P_V& predVec,
						std::map<std::string, double>& map_data_point ) const
{

  for( unsigned int i = 0; i < m_num_observables; ++i )
    {
      if( m_order_observables[i].obsLn )
	{
	  diffVec[i] = std::log( predVec[i] )
	    - std::log( map_data_point[ m_order_observables[i].obsName ] );
	}
      else
	{
	  diffVec[i] = predVec[i]
	    - map_data_point[ m_order_observables[i].obsName ];
	}
    }

}

//********************************************************************
// uqModel::get_correction_term
//********************************************************************
template <class P_V, class P_M>
  double uqModel<P_V,P_M>::get_correction_term( std::map<std::string, double>& map_data_point ) const
{

  double corr_term = 0.0;

  for( unsigned int i = 0; i < m_num_observables; ++i )
    {
      if( m_order_observables[i].obsLn )
	{
	  corr_term += std::log( map_data_point[ m_order_observables[i].obsName ] );
	}
    }

  return corr_term;

}

//********************************************************************
// uqModel::get_act_val_param
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::get_act_val_param( const P_V& paramValues,
					    P_V& actualValues )
{

  UQ_FATAL_TEST_MACRO( !m_set_params,
		       0,
		       "uqModel::get_act_val_param",
		       "!m_set_params - Computing actual values is not possible if all the parameters are not set" );
  
  for( unsigned int i = 0; i < m_num_params; ++i )
    {

      actualValues[ i ] = m_all_params[ m_order_params[ i ] ].getActualValue( paramValues[ i ] );

    }

}

//********************************************************************
// uqModel::test_printActualValue
//********************************************************************
template <class P_V, class P_M>
  void uqModel<P_V,P_M>::test_printActualValue( const P_V& paramValues )
{

  P_V actualValues( paramValues );

  get_act_val_param( paramValues, actualValues );

  std::cout << "Normalized param: ";
  paramValues.print( std::cout );
  std::cout << std::endl;
  std::cout << "Actual param: ";
  actualValues.print( std::cout );
  std::cout << std::endl;

}

#endif // __UQ_MODEL_H
