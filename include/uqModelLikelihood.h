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
 * uqModelLikelihood.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MODEL_LIKELIHOOD_H
#define __UQ_MODEL_LIKELIHOOD_H

#include "uqScalarFunction.h"
#include "uqGslMatrix.h"

template<class P_V,class P_M>
class uqModelLikelihood : public uqBaseScalarFunctionClass<P_V,P_M>
{

 public:
  
  uqModelLikelihood( const uqModel<P_V,P_M>& uq_model,
		     const uqDataSet& uq_data_set );

  virtual ~uqModelLikelihood();

  double actualValue( const P_V& domainVector, 
		      const P_V* domainDirection, 
		      P_V* gradVector, 
		      P_M* hessianMatrix, 
		      P_V* hessianEffect) const;

  double lnValue( const P_V& domainVector, 
		  const P_V* domainDirection, 
		  P_V* gradVector, 
		  P_M* hessianMatrix, 
		  P_V* hessianEffect) const;

  
 protected:

  uqModelLikelihood();

  const uqModel<P_V,P_M>& m_uq_model;
  const uqDataSet& m_uq_data_set;
  
};


//********************************************************************
// uqModelLikelihood::uqModelLikelihood
//********************************************************************
template <class P_V, class P_M>
  uqModelLikelihood<P_V,P_M>::uqModelLikelihood( const uqModel<P_V,P_M>& uq_model,
						 const uqDataSet& uq_data_set )
  : uqBaseScalarFunctionClass<P_V,P_M>( ("like_" + uq_model.getName()).c_str(), 
					*(uq_model.paramDomain) ),
  m_uq_model( uq_model ),
  m_uq_data_set( uq_data_set )
  {}

//********************************************************************
// uqModelLikelihood::~uqModelLikelihood
//********************************************************************
template <class P_V, class P_M>
  uqModelLikelihood<P_V,P_M>::~uqModelLikelihood()
{}

//********************************************************************
// uqModelLikelihood::actualValue
//********************************************************************
template<class P_V,class P_M>
  double uqModelLikelihood<P_V,P_M>::actualValue( const P_V& domainVector, 
						  const P_V* domainDirection, 
						  P_V* gradVector, 
						  P_M* hessianMatrix, 
						  P_V* hessianEffect) const
{

  double ln_value = m_uq_model.likelihoodFunction( domainVector,
						   domainDirection,
						   (void *) &m_uq_data_set,
						   gradVector,
						   hessianMatrix,
						   hessianEffect );

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  return std::exp(ln_value);
  
#else
  return std::exp(-.5*ln_value);
  
#endif

}

//********************************************************************
// uqModelLikelihood::lnValue
//********************************************************************
template<class P_V,class P_M>
  double uqModelLikelihood<P_V,P_M>::lnValue( const P_V& domainVector, 
					      const P_V* domainDirection, 
					      P_V* gradVector, 
					      P_M* hessianMatrix, 
					      P_V* hessianEffect) const
{

  double ln_value = this->m_uq_model.likelihoodFunction( domainVector,
							 domainDirection,
							 (void *) &(this->m_uq_data_set),
							 gradVector,
							 hessianMatrix,
							 hessianEffect );

  return ln_value;
}


#endif // __UQ_MODEL_LIKELIHOOD_H
