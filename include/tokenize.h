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
// tokenize.h: Utility function to tokenize a string and return a 
//             std::vector of the tokens.
//
// $Id: tokenize.h 12509 2010-08-18 22:12:34Z pbauman $ 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef TOKENIZE_H
#define TOKENIZE_H

#include <string>
#include <vector>

//********************************************************************
// Tokenize
//********************************************************************
void Tokenize(const std::string& str,
	      std::vector<std::string>& tokens,
	      const std::string& delimiters = " ")
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);

      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

#endif // TOKENIZE_H
