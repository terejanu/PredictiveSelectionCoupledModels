//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2008,2009 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of <APP/LIBRARY>.
//
// <APP/LIBRARY> is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// <APP/LIBRARY> is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with <APP/LIBRARY>.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
// 
// convert.h: Utility function for converting strings to doubles.
//
// $Id: convert.h 12509 2010-08-18 22:12:34Z pbauman $ 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef CONVERT_H
#define CONVERT_H

#include <string>

//********************************************************************
// StringToDouble
//********************************************************************
double StringToDouble(const std::string& string_in)
{
  std::istringstream is(string_in);

  double double_out;

  is >> double_out;

  return double_out;
}


#endif // CONVERT_H
