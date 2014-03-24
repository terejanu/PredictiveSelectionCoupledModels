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
 * uqBagModels.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_BAG_MODELS_H
#define __UQ_BAG_MODELS_H

#include <utility>
#include <fstream>
#include "uqVectorRV.h"
#include "uqInfoTheory.h"

// STRUCT: bagObject
template <class P_V, class P_M>
struct bagObject
{

  // model object
  uqModel<P_V,P_M>*                                     model;

  // inverse problem & postRV
  uqStatisticalInverseProblemClass<P_V,P_M>*            invProb;
  double                                                log_evidence;
  double                                                plausibility;
  uqGenericVectorRVClass<P_V,P_M>*                      postRV;

  // when the posterior is constructed from other sources
  uqBaseVectorSequenceClass<P_V,P_M>*                   postChain;
  uqSequentialVectorRealizerClass<P_V,P_M>*             postRealizer;

  // forward problem & qoiRV
  uqStatisticalForwardProblemClass< P_V,P_M,P_V,P_M >*  fwdProb;
  uqGenericVectorRVClass<P_V,P_M>*                      qoiRV;
  std::vector<double>                                   qoi_evidence;

  // obtain the exact samples used in the forward problem
  // equivalence needed to compute Mutual Information
  uqGenericVectorRVClass<P_V,P_M>*                      paramSmpRV;
  uqSequentialVectorRealizerClass<P_V,P_M>*             paramSmpRealizer;

};

// CLASS: uqBagModels
template <class P_V, class P_M>
class uqBagModels
{

 public:
  const uqBaseEnvironmentClass&                  m_env;
  std::vector< bagObject<P_V,P_M> >              m_bagModels;
  qoiDataStruct                                  m_qDS;
  uqVectorSpaceClass<P_V,P_M>*                   m_qoiSpace;
  uqGenericVectorRVClass<P_V,P_M>*               m_qoiRV;
  uqBaseVectorSequenceClass<P_V,P_M>*            m_qoiChain;
  uqSequentialVectorRealizerClass<P_V,P_M>*      m_qoiRealizer;

  uqBagModels( const uqBaseEnvironmentClass& env );
  uqBagModels( const uqBaseEnvironmentClass& env,
	       std::string bagName );
  ~uqBagModels();

  void addModel( uqModel<P_V,P_M>* my_model );
  void addQoI( const qoiDataStruct& qDS );
  void compute_model_plausibility();
  void compute_combined_qoi();
  void compute_qoi_aware_ev();
  void save_combined_qoi();
  void compute_conditional_mi( std::string true_qoi, 
			       std::vector<std::string>& surrogate_qois );
  void print_evidences();
  void print_conditional_mi();

 private:
  std::string                                    m_prefix;
  bool                                           m_qoi_set;
  std::string                                    m_true_qoi;
  std::vector<std::string>                       m_surrogate_qois;
  std::vector<double>                            m_surrogate_qois_cond_mi;

  void initialize_prediction_var();

};


//********************************************************************
// uqBagModels::uqBagModels
//********************************************************************
template <class P_V, class P_M>
  uqBagModels<P_V,P_M>::uqBagModels( const uqBaseEnvironmentClass& env )
  : m_prefix( "NONE" ),
  m_env( env ),
  m_qoi_set( false ),
  m_qoiSpace( NULL ),
  m_qoiChain( NULL),
  m_qoiRealizer( NULL ),
  m_qoiRV( NULL )
  {}

//********************************************************************
// uqBagModels::uqBagModels
//********************************************************************
template <class P_V, class P_M>
  uqBagModels<P_V,P_M>::uqBagModels( const uqBaseEnvironmentClass& env,
				     std::string bagName )
  : m_prefix( bagName ),
  m_env( env ),
  m_qoi_set( false ),
  m_qoiSpace( NULL ),
  m_qoiChain( NULL),
  m_qoiRealizer( NULL ),
  m_qoiRV( NULL )
  {}

//********************************************************************
// uqBagModels::~uqBagModels
//********************************************************************
template <class P_V, class P_M>
  uqBagModels<P_V,P_M>::~uqBagModels()
{

  for( unsigned int i = 0; i < m_bagModels.size(); ++i )
    {
      //      std::cout << m_prefix << " 1" << std::endl;
      if( m_bagModels[i].paramSmpRealizer )   delete m_bagModels[i].paramSmpRealizer;
      //      std::cout << m_prefix << " 2" << std::endl;
      if( m_bagModels[i].postChain ) 	      delete m_bagModels[i].postChain;
      //      std::cout << m_prefix << " 2" << std::endl;
      if( m_bagModels[i].postRealizer )	      delete m_bagModels[i].postRealizer;
      //      std::cout << m_prefix << " 2" << std::endl;
      if( m_bagModels[i].postRV ) 	      delete m_bagModels[i].postRV;
      //      std::cout << m_prefix << " 3" << std::endl;
      if( m_bagModels[i].qoiRV ) 	      delete m_bagModels[i].qoiRV;
      //      std::cout << m_prefix << " 4" << std::endl;
      if( m_bagModels[i].paramSmpRV )         delete m_bagModels[i].paramSmpRV;
      //      std::cout << m_prefix << " 5" << std::endl;
      if( m_bagModels[i].fwdProb )            delete m_bagModels[i].fwdProb;
      //      std::cout << m_prefix << " 6" << std::endl;
      if( m_bagModels[i].invProb )            delete m_bagModels[i].invProb;
    }

  //  std::cout << m_prefix << " 7" << std::endl;
  if( m_qoiChain )                            delete m_qoiChain;
  //  std::cout << m_prefix << " 8" << std::endl;
  if( m_qoiRealizer )                         delete m_qoiRealizer;
  //  std::cout << m_prefix << " 9" << std::endl;
  if( m_qoiSpace )                            delete m_qoiSpace;
  //  std::cout << m_prefix << " 10" << std::endl;
  if( m_qoiRV )                               delete m_qoiRV;
  //  std::cout << m_prefix << " 11" << std::endl;

}

//********************************************************************
// uqBagModels::addModel
//********************************************************************
template <class P_V, class P_M>
  void uqBagModels<P_V,P_M>::addModel( uqModel<P_V,P_M>* my_model )
{

  bagObject<P_V,P_M> newModelObject;

  newModelObject.model = my_model;
  newModelObject.postRV = new uqGenericVectorRVClass<P_V,P_M>
    ( ("post_" + my_model->getName()).c_str(), *(my_model->paramSpace) );
  newModelObject.invProb = NULL;
  newModelObject.fwdProb = NULL;
  newModelObject.postChain = NULL;
  newModelObject.postRealizer = NULL;

  if( m_qoi_set )
    {
      newModelObject.qoiRV = new uqGenericVectorRVClass<P_V,P_M>
	    ( ("qoi_" + my_model->getName()).c_str(), *m_qoiSpace );

      newModelObject.paramSmpRV = new uqGenericVectorRVClass<P_V,P_M>
	( ("postSmp_" + my_model->getName()).c_str(), *(my_model->paramSpace) );  

      newModelObject.paramSmpRealizer = NULL;
    }
  else
    {
      newModelObject.qoiRV = NULL;
      newModelObject.paramSmpRV = NULL;
      newModelObject.paramSmpRealizer = NULL;
    }

  m_bagModels.push_back( newModelObject );  
  
}

//********************************************************************
// uqBagModels::addQoI
//********************************************************************
template <class P_V, class P_M>
  void uqBagModels<P_V,P_M>::addQoI( const qoiDataStruct& qDS )
{

  m_qDS = qDS;
  m_qoi_set = true;

  m_qoiSpace = new uqVectorSpaceClass<P_V,P_M>( m_env, "qoi_", m_qDS.m_qois.size(), NULL );

  // initialize the qoiRV for all the models so far
  if( m_bagModels.size() > 0 )
    {
      initialize_prediction_var();
    }

}

//********************************************************************
// uqBagModels::initialize_prediction_var
//********************************************************************
template <class P_V, class P_M>
  void uqBagModels<P_V,P_M>::initialize_prediction_var( )
{

  for( unsigned int i = 0; i < m_bagModels.size(); ++i )
    {
      if( !m_bagModels[i].qoiRV )
	{

	  m_bagModels[i].qoiRV = new uqGenericVectorRVClass<P_V,P_M>
	    ( ("qoi_" + m_bagModels[i].model->getName()).c_str(), *m_qoiSpace );

	  m_bagModels[i].paramSmpRV = new uqGenericVectorRVClass<P_V,P_M>
	    ( ("paramSmp_" + m_bagModels[i].model->getName()).c_str(), *(m_bagModels[i].model->paramSpace) );  

	}
    } 

}

//********************************************************************
// uqBagModels::compute_model_plausibility
//********************************************************************
template <class P_V, class P_M>
  void  uqBagModels<P_V,P_M>::compute_model_plausibility()
{

  double total_pl = 0.0;

  for( unsigned int i = 0; i < m_bagModels.size(); ++i )
    {
      total_pl += std::exp( m_bagModels[i].log_evidence );
    }  

  for( unsigned int i = 0; i < m_bagModels.size(); ++i )
    {
      m_bagModels[i].plausibility = std::exp( m_bagModels[i].log_evidence ) / total_pl;
    }

}


//********************************************************************
// uqBagModels::compute_combined_qoi
// (the model plausibilities has to be computed first)
//********************************************************************
template <class P_V, class P_M>
  void  uqBagModels<P_V,P_M>::compute_combined_qoi()
{

  // get the min size qoi chain in the bag
  unsigned int sizeQoI = m_bagModels[0].qoiRV->realizer().subPeriod();

  for( unsigned int i = 0; i < m_bagModels.size(); ++i )
    {
      if( sizeQoI > m_bagModels[i].qoiRV->realizer().subPeriod() )
	{
	  sizeQoI = m_bagModels[i].qoiRV->realizer().subPeriod();
	}
    }

  std::cout << "Combined QoI will have " << sizeQoI << " samples" << std::endl;

  // generate the m_qoiChain
  m_qoiChain = new uqSequenceOfVectorsClass<P_V,P_M>( *m_qoiSpace, sizeQoI, "qoiChain_"+m_prefix );
  P_V smpVec( m_qoiSpace->zeroVector() );

  unsigned int lastIdx = 0;
  for( unsigned int i = 0; i < m_bagModels.size(); ++i )
    {

      // get the number we want to use from each chain by weighting
      // it using plausibility. Two variants here:
      //    - The new RV has the same number for samples as the rest of the samples
      //    - The new RV has more samples than the rest
      unsigned int i_sizeQoI;
      if( (i+1) == m_bagModels.size() )
	{
	  i_sizeQoI = sizeQoI - lastIdx;
	}
      else
	{
	  i_sizeQoI = (unsigned int) floor( (double)m_bagModels[i].qoiRV->realizer().subPeriod() * m_bagModels[i].plausibility );
	}

      std::cout << "Model " << i << " qoi smp = " << m_bagModels[i].qoiRV->realizer().subPeriod() << std::endl;
      std::cout << "Model " << i << " plausibility = " << m_bagModels[i].plausibility << std::endl;
      std::cout << "Model " << i << " has " << i_sizeQoI << " samples to combine" << std::endl;

      for( unsigned int j = 0; j < i_sizeQoI; ++j )
	{
	  m_bagModels[i].qoiRV->realizer().realization( smpVec );
	  m_qoiChain->setPositionValues( lastIdx + j, smpVec );
	}
      lastIdx += i_sizeQoI;

    }

  // instantiate the m_qoiRV
  m_qoiRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>
    ( ("qoiRealizer_"+m_prefix).c_str(), *m_qoiChain );

  m_qoiRV = new uqGenericVectorRVClass<P_V,P_M>( ("qoi_"+m_prefix).c_str(), *m_qoiSpace );
  m_qoiRV->setRealizer( *m_qoiRealizer );

}

//********************************************************************
// uqBagModels::print_evidences
//********************************************************************
template <class P_V, class P_M>
  void  uqBagModels<P_V,P_M>::print_evidences()
{

  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << std::setw(20) << "Model Name" << 
    std::setw(20) << "Log Evidence" << 
    std::setw(20) << "Plausibility";
  for( unsigned int i = 0; i < m_qDS.m_qois.size(); ++i )
    {
      std::cout << std::setw(20) << m_qDS.m_qois[i].substr(0,16);
    }
  std::cout << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;

  for( unsigned int j = 0; j < m_bagModels.size(); ++j )
    {

      std::cout << std::scientific << std::setw(20) <<
	m_bagModels[j].model->getName().substr(0,16) << 
	std::setw(20) << m_bagModels[j].log_evidence << 
	std::setw(20) << m_bagModels[j].plausibility;
      for( unsigned int i = 0; i < m_bagModels[j].qoi_evidence.size(); ++i )
	{
	  std::cout << std::scientific << std::setw(20) << 
	    m_bagModels[j].qoi_evidence[i];
	}
      std::cout << std::endl;

    }

  std::cout << "-----------------------------------------------------" << std::endl;

}

//********************************************************************
// uqBagModels::compute_qoi_aware_ev
// NOTE: all the inverse and forward problems have to be ran first
//       also the combined QoI has to be estimated
//********************************************************************
template <class P_V, class P_M>
  void uqBagModels<P_V,P_M>::compute_qoi_aware_ev()
{

  unsigned int k = 1;
  double eps = UQ_INFTH_ANN_EPS;
  unsigned int xN = m_qoiRV->realizer().subPeriod();
  unsigned int dimX = 1;
  unsigned int dimY = 1;
  unsigned int xDimSel[1] = { 0 };
  unsigned int yDimSel[1] = { 0 }; 

  for( unsigned int idx = 0; idx < m_bagModels.size(); ++idx )
    {
      unsigned int yN = m_bagModels[idx].qoiRV->realizer().subPeriod();

      for( unsigned int i = 0; i < m_qDS.m_qois.size(); ++i )
	{
	  
	  xDimSel[0] = i;
	  yDimSel[0] = i;
	  
	  /*
	    double tmp = -estimateCE_ANN( *m_qoiRV, 
	    *(m_bagModels[idx].qoiRV), 
	    xDimSel, dimX,
	    yDimSel, dimY,
	    xN, yN,
	    k, eps );
	  */
	  
	  double tmp = estimateKL_ANN( *m_qoiRV, 
				       *(m_bagModels[idx].qoiRV), 
				       xDimSel, dimX,
				       yDimSel, dimY,
				       xN, yN,
				       k, eps );

	  m_bagModels[idx].qoi_evidence.push_back( tmp );

	}
    }

}


//********************************************************************
// uqBagModels::save_combined_qoi
//********************************************************************
template <class P_V, class P_M>
  void  uqBagModels<P_V,P_M>::save_combined_qoi()
{

  std::string fileName = "./outputData/" + m_prefix + "_qoi.m";
  unsigned int rowsQoI = m_qoiRV->realizer().subPeriod();
  unsigned int colsQoI = m_qDS.m_qois.size();
  P_V smpVec( m_qoiSpace->zeroVector() );

  std::ofstream qoiFile;
  qoiFile.open ( fileName.c_str() );

  qoiFile << "combined_qoi_seq = zeros(" << rowsQoI << "," << colsQoI << ");" << std::endl;
  qoiFile << "combined_qoi_seq = [";

  for( unsigned int i = 0; i < rowsQoI; ++i )
    {
      m_qoiRV->realizer().realization( smpVec );
      for( unsigned int j = 0; j < colsQoI-1; ++j )
	{
	  qoiFile << std::scientific << smpVec[j] << " ";
	}
      qoiFile << std::scientific << smpVec[colsQoI-1] << std::endl;
    }
  
  qoiFile << "];" << std::endl;

  qoiFile.close();

}

//********************************************************************
// uqBagModels::compute_conditional_mi
//********************************************************************
template <class P_V, class P_M>
  void  uqBagModels<P_V,P_M>::compute_conditional_mi( std::string true_qoi, 
						      std::vector<std::string>& surrogate_qois )
{

  // set member variables
  m_true_qoi = true_qoi;
  m_surrogate_qois = surrogate_qois;

  // prepare calculations
  m_surrogate_qois_cond_mi.clear();

  unsigned int dimX = 1;                // dim of the true qoi
  unsigned int dimY = 1;                // dim of each surrogate qoi

  unsigned int xDimSel[1] = { 0 };      // position true qoi in the joint
  unsigned int yDimSel[1] = { 0 };      // position surrogate qoi in the joint
  std::vector<std::string>::iterator i; 

  unsigned int k = UQ_INFTH_ANN_KNN;                 // k nearest neighbor
  unsigned int eps = UQ_INFTH_ANN_EPS;               // tolerance in finding kNN

  // find position of the true qoi
  i = find(m_qDS.m_qois.begin(), m_qDS.m_qois.end(), m_true_qoi);

  UQ_FATAL_TEST_MACRO( i == m_qDS.m_qois.end(),
		       0,
		       "uqBagModels::compute_conditional_mi",
		       "could not find the true qoi" );

  xDimSel[0] = i - m_qDS.m_qois.begin();

  // compute conditional mutual information for all tst qois
  for( unsigned int i_qoi = 0; i_qoi < m_surrogate_qois.size(); ++i_qoi )
    {

      double tmp_cond_mi = 0.0;

      // find position of the surrogate qoi
      i = find(m_qDS.m_qois.begin(), m_qDS.m_qois.end(), m_surrogate_qois[i_qoi]);

      UQ_FATAL_TEST_MACRO( i == m_qDS.m_qois.end(),
			   0,
			   "uqBagModels::compute_conditional_mi",
			   "could not find the surrogate qoi" );

      yDimSel[0] = i - m_qDS.m_qois.begin();

      // loop over all the models
      for( unsigned int j_model = 0; j_model < m_bagModels.size(); ++j_model )
	{

	  // get the number of samples
	  unsigned int N = m_bagModels[j_model].qoiRV->realizer().subPeriod();    

	  double tmp_MI = estimateMI_ANN( *(m_bagModels[j_model].qoiRV),
					  xDimSel, dimX,
					  yDimSel, dimY,
					  k, N, eps ); 

	  std::cout << m_bagModels[j_model].model->getName() << " : MI<" <<
	    m_surrogate_qois[i_qoi] << "> = " << tmp_MI << std::endl;

	  tmp_cond_mi += tmp_MI * m_bagModels[j_model].plausibility;

	}

      // save conditional mutual information
      m_surrogate_qois_cond_mi.push_back( tmp_cond_mi );

    }


}

//********************************************************************
// uqBagModels::print_conditional_mi
//********************************************************************
template <class P_V, class P_M>
  void  uqBagModels<P_V,P_M>::print_conditional_mi()
{

  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "CONDITIONAL MUTUAL INFORMATION " << std::endl;
  std::cout << "Name Bag = " << m_prefix << std::endl;
  for( unsigned int i = 0; i < m_surrogate_qois.size(); ++i )
    {
      std::cout << "I( " << m_true_qoi << " ; " << m_surrogate_qois[i] << " ) = " <<
	m_surrogate_qois_cond_mi[i] << std::endl;
    }
  std::cout << "-----------------------------------------------------" << std::endl;

}


#endif // __UQ_BAG_MODELS_H


