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
 * uqParameter.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_PARAMETER_H
#define __UQ_PARAMETER_H

#include <math.h>
#include "uqEnvironment.h"

class uqParameter
{
 public:

  uqParameter();
  uqParameter( std::string name );
  uqParameter( std::string name, 
	       double min_val,
	       double max_val,
	       bool log_scale );
  ~uqParameter();

  void setName( std::string name );
  void setParam( double min_val,
		 double max_val,
		 bool log_scale );
  void getParam( double &min_val,
		 double &max_val,
		 bool &log_scale ) const;

  void setParam( std::string name, 
		 double min_val,
		 double max_val,
		 bool log_scale );

  // TODO: add a nominal value / normalized nominal value

  double getNormMin();
  double getNormMax();
  double getActualValue( double norm_value ) const;

  void printInfo();

 private:

  std::string _name;

  double _min;
  double _max;
  double _norm_val;
  double _norm_min;
  double _norm_max;
  bool _log_scale;

  void normalizeParam();

};

//********************************************************************
// uqParameter::uqParameter
//********************************************************************
uqParameter::uqParameter()
: _name( "NONE" )
{}

//********************************************************************
// uqParameter::uqParameter
//********************************************************************
uqParameter::uqParameter( std::string name )
: _name( name )
{}
  
//********************************************************************
// uqParameter::uqParameter
//********************************************************************
uqParameter::uqParameter( std::string name,
			  double min_val,
			  double max_val,
			  bool log_scale )
: _name( name ),
  _min( min_val ),
  _max( max_val ),
  _log_scale( log_scale )
{

  normalizeParam();
  
}

//********************************************************************
// uqParameter::~uqParameter
//********************************************************************
uqParameter::~uqParameter()
{}

//********************************************************************
// uqParameter::setParam
//********************************************************************
void uqParameter::setName( std::string name )
{

  _name = name;
  
}

//********************************************************************
// uqParameter::setParam
//********************************************************************
void uqParameter::setParam( double min_val,
			    double max_val,
			    bool log_scale )
{

  _min = min_val;
  _max = max_val;
  _log_scale = log_scale;

  normalizeParam();
  
}

//********************************************************************
// uqParameter::getParam
//********************************************************************
void uqParameter::getParam( double &min_val,
			    double &max_val,
			    bool &log_scale ) const
{

  min_val = _min;
  max_val = _max;
  log_scale = _log_scale;

}


//********************************************************************
// uqParameter::setParam
//********************************************************************
void uqParameter::setParam( std::string name, 
			    double min_val,
			    double max_val,
			    bool log_scale )
{

  _name = name;
  _min = min_val;
  _max = max_val;
  _log_scale = log_scale;

  normalizeParam();
  
}

//********************************************************************
// uqParameter::getNormMin
//********************************************************************
double uqParameter::getNormMin( )
{
  return _norm_min;
}

//********************************************************************
// uqParameter::getNormMax
//********************************************************************
double uqParameter::getNormMax( )
{
  return _norm_max;
}

//********************************************************************
// uqParameter::getActualValue
//********************************************************************
double uqParameter::getActualValue( double norm_value ) const
{

  double actual_value;

  actual_value = norm_value * _norm_val;
  
  if( _log_scale )
    {
      actual_value = pow( 10.0, actual_value );
    }

  return actual_value;

}

//********************************************************************
// uqParameter::printInfo
//********************************************************************
void uqParameter::printInfo()
{

  std::cout << "Param: " << _name << std::endl;
  std::cout << "\t Minim     = " << _min << std::endl;
  std::cout << "\t Maxim     = " << _max << std::endl;
  std::cout << "\t Log.Scale = " << _log_scale << std::endl;
  std::cout << "\t Norm.Val  = " << _norm_val << std::endl;
  std::cout << "\t Norm.Min  = " << _norm_min << std::endl;
  std::cout << "\t Norm.Max  = " << _norm_max << std::endl;

}

//********************************************************************
// uqParameter::normalizeParam
//********************************************************************
void uqParameter::normalizeParam()
{

  UQ_FATAL_TEST_MACRO( _max < _min,
		       0,
		       "uqParameter::normalizeParam",
		       "_max < _min" );

  _norm_val = _max - _min;

  UQ_FATAL_TEST_MACRO( _norm_val == 0,
		       0,
		       "uqParameter::normalizeParam",
		       "_norm_val == 0" );

  _norm_min = _min / _norm_val;
  _norm_max = _max / _norm_val;
  
}

#endif // __UQ_PARAMETER_H
