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
 * uqStatsAlg.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_STATS_ALG_H
#define __UQ_STATS_ALG_H

#include "stdio.h"
#include "uqModel.h"
#include "uqDataSet.h"
#include "uqModelLikelihood.h"
#include "uqModelQoI.h"
#include "uqStatisticalInverseProblem.h"
#include "uqStatisticalForwardProblem.h"
#include "uqBagModels.h"
#include "uqInfoTheory.h"

//********************************************************************
// FUNC.DEF.
//********************************************************************
void selectiveRename(std::string s_file, std::string m_prefix);

template <class P_V, class P_M>
  void renameFilesIP( uqModel<P_V,P_M>& my_model );

//********************************************************************
// solveIP
//********************************************************************
template <class P_V, class P_M>
  double solveIP( bagObject<P_V,P_M>& my_model_obj, uqDataSet& my_data_set )
{

  std::string prefix_ip = "";

  // define the model likelihood
  uqModelLikelihood<P_V,P_M> myLikelihood( *(my_model_obj.model), my_data_set );

  // defne a new statistical inverse problem
  my_model_obj.invProb = new uqStatisticalInverseProblemClass<P_V,P_M>
    ( prefix_ip.c_str(), 
      NULL,
      *(my_model_obj.model->priorRV),
      myLikelihood, 
      *(my_model_obj.postRV) );

  // solve the statistical inverse problem
  my_model_obj.invProb->solveWithBayesMLSampling();

  // print out
  std::cout << "********* INFO INVERSE PROBLEM ***********" << std::endl;
  std::cout << "nameModel   = " << my_model_obj.model->getName() << std::endl;
  std::cout << "logEvidence = " << my_model_obj.invProb->logEvidence() << std::endl;
  std::cout << "evidence    = " << exp(my_model_obj.invProb->logEvidence()) << std::endl;
  std::cout << "expLogData  = " << my_model_obj.invProb->meanLogLikelihood() << std::endl;
  std::cout << "infoGain    = " << my_model_obj.invProb->eig() << std::endl;

  // try to compute this infoGain from Information Theoretic measures
  // assume that the prior is uniform and normalized (volume = 1) as in this case
  double estInfoGain = - my_model_obj.postRV->estimateENT_ANN();
  std::cout << "estInfoGain = " << estInfoGain << std::endl;
  std::cout << "******************************************" << std::endl;

  // save evidence
  my_model_obj.log_evidence = my_model_obj.invProb->logEvidence();

  // rename the files generated
  renameFilesIP( *(my_model_obj.model) );

  return my_model_obj.invProb->logEvidence();

}

//********************************************************************
// solveFP
// NOTE: this uses postRV to get the QoIs, so the inverse problem
//       has to be solved first
//********************************************************************
template <class P_V, class P_M>
  void solveFP( uqBagModels<P_V,P_M>& allModels, unsigned int i )
{

  std::string prefix_fp = "";

  bagObject<P_V,P_M>* my_model_qoi = &(allModels.m_bagModels[i]);

  // define the model QoI
  uqModelQoI< P_V,P_M,P_V,P_M > 
    myQoI( *(my_model_qoi->model), 
	   allModels.m_qDS, 
	   my_model_qoi->qoiRV->imageSet() );

  // define a new statistical forward problem
  my_model_qoi->fwdProb = new uqStatisticalForwardProblemClass< P_V,P_M,P_V,P_M > 
    ( prefix_fp.c_str(), 
      NULL,
      *(my_model_qoi->postRV), 
      myQoI, 
      *(my_model_qoi->qoiRV) );

  // solve the forward problem
  my_model_qoi->fwdProb->solveWithMonteCarlo( NULL );

  // set the realizer to obtain the exact param samples that generated the QoI samples
  my_model_qoi->paramSmpRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>
    (("paramRealizer_" + my_model_qoi->model->getName()).c_str(), 
     my_model_qoi->fwdProb->getParamChain() );

  my_model_qoi->paramSmpRV->setRealizer( *(my_model_qoi->paramSmpRealizer) );

  // rename the files generated
  renameFilesFP( *(my_model_qoi->model) );

}



//********************************************************************
// constructCoupledPostRV
// NOTE: this uses postRV from component models, so these inverse 
//       problems have to be solved first
//********************************************************************
template <class P_V, class P_M>
  void constructCoupledPostRV( uqBagModels<P_V,P_M>& coupledModels, 
			       std::string coupledName,
			       uqBagModels<P_V,P_M>& principalModels, 
			       std::string principalName,
			       uqBagModels<P_V,P_M>& secondaryModels, 
			       std::string secondaryName )
{

  // find for the coupled model in the bag
  bagObject<P_V,P_M>* my_coupled_object;
  for( unsigned int i = 0; i < coupledModels.m_bagModels.size(); ++i )
    {
      if( coupledModels.m_bagModels[i].model->getName() == coupledName )
	{
	  my_coupled_object = &(coupledModels.m_bagModels[i]);
	  break;
	}
    }

  // find for the principal model in the bag
  bagObject<P_V,P_M>* my_principal_object;
  for( unsigned int i = 0; i < principalModels.m_bagModels.size(); ++i )
    {
      if( principalModels.m_bagModels[i].model->getName() == principalName )
	{
	  my_principal_object = &(principalModels.m_bagModels[i]);
	  break;
	}
    }

  // find for the secondary model in the bag
  bagObject<P_V,P_M>* my_secondary_object;
  for( unsigned int i = 0; i < secondaryModels.m_bagModels.size(); ++i )
    {
      if( secondaryModels.m_bagModels[i].model->getName() == secondaryName )
	{
	  my_secondary_object = &(secondaryModels.m_bagModels[i]);
	  break;
	}
    }

  // sanity check: dim(coupled postRV) = dim(principal postRV) + dim(secondary postRV)
  unsigned int coupled_dim = my_coupled_object->model->getNumberParams();
  unsigned int principal_dim = my_principal_object->model->getNumberParams();
  unsigned int secondary_dim = my_secondary_object->model->getNumberParams();

  UQ_FATAL_TEST_MACRO( coupled_dim != (principal_dim + secondary_dim),
		       0,
		       "uqStatsAlg.h : constructCoupledPostRV",
		       "dim(coupled postRV) != dim(principal postRV) + dim(secondary postRV)" );
  
  // sanity check: the size of the two chains should also agree and 
  // it should greater than zero
  // NOTE: that this can be relaxed
  unsigned int principal_chain_size = my_principal_object->postRV->realizer().subPeriod();
  unsigned int secondary_chain_size = my_secondary_object->postRV->realizer().subPeriod();

  UQ_FATAL_TEST_MACRO( principal_chain_size != secondary_chain_size,
		       0,
		       "uqStatsAlg.h : constructCoupledPostRV",
		       "the size of the two chain does not match" );

  UQ_FATAL_TEST_MACRO( principal_chain_size == 0,
		       0,
		       "uqStatsAlg.h : constructCoupledPostRV",
		       "chain size is zero, which indicates that the Inverse Problems have not been run" );

  unsigned int coupled_chain_size = principal_chain_size;

  // create the three vectors for retriving the chain samples
  P_V coupledSmpVec( my_coupled_object->model->paramSpace->zeroVector() );
  P_V principalSmpVec( my_principal_object->model->paramSpace->zeroVector() );
  P_V secondarySmpVec( my_secondary_object->model->paramSpace->zeroVector() );

  // instantiate chain
  my_coupled_object->postChain = new uqSequenceOfVectorsClass<P_V,P_M>
    ( my_coupled_object->postRV->imageSet().vectorSpace(),
      coupled_chain_size,
      my_coupled_object->model->getName() + "paramChain" );
  
  // populate chain
  for( unsigned int i = 0; i < coupled_chain_size; ++i )
    {

      // get samples
      my_principal_object->postRV->realizer().realization( principalSmpVec );
      my_secondary_object->postRV->realizer().realization( secondarySmpVec );
      
      // concatenate vectors and add to the chain
      coupledSmpVec.cwSetConcatenated( principalSmpVec, secondarySmpVec );
      my_coupled_object->postChain->setPositionValues( i, coupledSmpVec );
      
    }

  // instantiate and set realizer
  my_coupled_object->postRealizer = new uqSequentialVectorRealizerClass<P_V,P_M>
    (("paramRealizer_" + my_coupled_object->model->getName()).c_str(), 
     *(my_coupled_object->postChain) );
  
  // attach realizer to the postRV
  my_coupled_object->postRV->setRealizer( *(my_coupled_object->postRealizer) );

  // compute and save Model Plausibility for the coupled model
  // p(M^AB|D^A,D^B) = P(M^A|D^A) * P(M^B|D^B)
  my_coupled_object->plausibility = my_principal_object->plausibility *
    my_secondary_object->plausibility;

  // compute and save the evidence for the coupled model
  // log( p(D^A,D^B|M^AB) ) = log( P(D^A|M^A) ) + log( P(D^B|M^B) )
  my_coupled_object->log_evidence = my_principal_object->log_evidence +
    my_secondary_object->log_evidence;

}

//********************************************************************
// renameFiles : IP
//********************************************************************
template <class P_V, class P_M>
  void renameFilesIP( uqModel<P_V,P_M>& my_model ) 
{
  // TODO: these files names should be obtained from queso.inp
  selectiveRename( "rawChain_ml.m", my_model.getName() );
  selectiveRename( "sipOutput_ml_sub0.m", my_model.getName() );
}

//********************************************************************
// renameFiles : FP
//********************************************************************
template <class P_V, class P_M>
  void renameFilesFP( uqModel<P_V,P_M>& my_model ) 
{
  // TODO: these files names should be obtained from queso.inp
  selectiveRename( "fp_p_seq.m", my_model.getName() );
  selectiveRename( "fp_p_seq_sub0.m", my_model.getName() );
  selectiveRename( "fp_q_seq.m", my_model.getName() );
  selectiveRename( "fp_q_seq_sub0.m", my_model.getName() );
  selectiveRename( "sfpOutput_sub0.m", my_model.getName() );
}

//********************************************************************
// selectiveRename
//********************************************************************
void selectiveRename(std::string s_file, std::string m_prefix) 
{

  std::string file_s1 = "./outputData/" + s_file;
  std::string file_s2 = "./outputData/" + m_prefix + "_" + s_file;

  int resTmp = rename( file_s1.c_str(), file_s2.c_str() );

  UQ_FATAL_TEST_MACRO( resTmp != 0,
		       0,
		       "uqStatsAlg - selectiveRename",
		       "Could not rename files" );
  
}



#endif // __UQ_STATS_ALG_H

