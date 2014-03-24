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
 * uqModelQoI.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_MODEL_QOI_H
#define __UQ_MODEL_QOI_H

#include "uqVectorFunction.h"
#include "uqGslMatrix.h"
#include "uqDistArray.h"
#include "uqModel.h"

template<class P_V,class P_M,class Q_V,class Q_M>
  class uqModelQoI : public uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>
{

 public:
  
  uqModelQoI( const uqModel<P_V,P_M>& uq_model,
	      const qoiDataStruct& qDS,
	      const uqVectorSetClass<Q_V,Q_M>& qoi_imageSet );

  virtual ~uqModelQoI();

  void compute ( const P_V&              domainVector,
		 const P_V*              domainDirection,
		 Q_V&                    imageVector,
		 uqDistArrayClass<P_V*>* gradVectors,     // Yes, 'P_V'
		 uqDistArrayClass<P_M*>* hessianMatrices, // Yes, 'P_M'
		 uqDistArrayClass<P_V*>* hessianEffects ) const;

 protected:

  uqModelQoI();

  const uqModel<P_V,P_M>& m_uq_model;
  const qoiDataStruct& m_qDS;

};


//********************************************************************
// uqModelQoI::uqModelQoI
//********************************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
  uqModelQoI<P_V,P_M,Q_V,Q_M>::uqModelQoI( const uqModel<P_V,P_M>& uq_model,
					   const qoiDataStruct& qDS,
					   const uqVectorSetClass<Q_V,Q_M>& qoi_imageSet )
  : uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>( ("qoi_" + uq_model.getName()).c_str(), 
						uq_model.priorRV->imageSet(),
						qoi_imageSet ),
  m_uq_model( uq_model ),
  m_qDS( qDS )
{}

//********************************************************************
// uqModelQoI::~uqModelQoI
//********************************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
  uqModelQoI<P_V,P_M,Q_V,Q_M>::~uqModelQoI()
{}

//********************************************************************
// uqModelQoI::compute
//********************************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
  void uqModelQoI<P_V,P_M,Q_V,Q_M>::compute( const P_V&              domainVector,
					     const P_V*              domainDirection,
					     Q_V&                    imageVector,
					     uqDistArrayClass<P_V*>* gradVectors,     // Yes, 'P_V'
					     uqDistArrayClass<P_M*>* hessianMatrices, // Yes, 'P_M'
					     uqDistArrayClass<P_V*>* hessianEffects ) const
{

  m_uq_model.qoiFunction( domainVector,
			  domainDirection,
			  (void *) &m_qDS,
			  imageVector,
			  gradVectors,
			  hessianMatrices,
			  hessianEffects );

}

#endif // __UQ_MODEL_QOI_H
