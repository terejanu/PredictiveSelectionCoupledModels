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
 * uqDataSet.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_DATASET_H
#define __UQ_DATASET_H

#include <map>
#include "uqDefines.h"
#include "tokenize.h"
#include "convert.h"

class uqDataSet
{

 public:

  uqDataSet();
  uqDataSet( std::string fileName );
  ~uqDataSet();

  unsigned int getNoColumns();
  unsigned int getDataSetSize() const;
  std::string getFileName();

  void readFromFile( std::string fileName );
  double getValue( std::string nameVar, unsigned int idx );
  void getMapDataPoint( unsigned int idx, std::map<std::string,double>& map_data_point ) const;

  void printInfo() const;

 private:

  unsigned int                                   m_data_set_size;
  unsigned int                                   m_no_columns;
  std::string                                    m_file_name;

  std::map< std::string,std::vector<double> >    m_data_column;
  std::vector<std::string>                       m_column_name;

  void cleanContent();

};


//********************************************************************
// uqDataSet::uqDataSet
//********************************************************************
uqDataSet::uqDataSet()
: m_data_set_size( 0 ),
  m_no_columns( 0 ),
  m_file_name( "" )
{}

//********************************************************************
// uqDataSet::uqDataSet
//********************************************************************
uqDataSet::uqDataSet( std::string fileName )
: m_data_set_size( 0 ),
  m_no_columns( 0 ),
  m_file_name( "" )
{
  readFromFile( fileName );
}

//********************************************************************
// uqDataSet::~uqDataSet
//********************************************************************
uqDataSet::~uqDataSet()
{}

//********************************************************************
// uqDataSet::getNoColumns
//********************************************************************
unsigned int uqDataSet::getNoColumns()
{

  return m_no_columns;

}

//********************************************************************
// uqDataSet::getDataSetSize
//********************************************************************
unsigned int uqDataSet::getDataSetSize() const
{

  return m_data_set_size;

}

//********************************************************************
// uqDataSet::getFileName
//********************************************************************
std::string uqDataSet::getFileName()
{

  return m_file_name;

}

//********************************************************************
// uqDataSet::readFromFile
//********************************************************************
void uqDataSet::readFromFile( std::string fileName )
{

  cleanContent();

  std::ifstream data_in( fileName.c_str(), std::ios::in );

  UQ_FATAL_TEST_MACRO( !data_in.is_open(),
		       0,
		       "uqDataSet::readFromFile",
		       "Could not open file for reading" );
  
  std::string row;
  std::vector<std::string> tokens;

  // Read first line of file to obtain the name of the columns
  getline(data_in, row);

  UQ_FATAL_TEST_MACRO( row.empty(),
		       0,
		       "uqDataSet::readFromFile",
		       "The header line is empty - no column names" );

  // Be sure to clear tokens each time since Tokenize uses push_back
  tokens.clear();
  Tokenize( row, tokens, " " );
  
  for( unsigned int i = 0; i < tokens.size(); ++i )
    {      
      m_data_column[ tokens[i] ].clear();
      m_column_name.push_back( tokens[i] );
    }

  m_no_columns = m_column_name.size();

  // Now read the rest of the data until we reach the end of the file
  while( !data_in.eof() )
    {
      getline(data_in, row);

      // Be sure to clear tokens each time since Tokenize uses push_back
      tokens.clear();

      if( !row.empty() ) 
	{
	  Tokenize( row, tokens, " " );

	  UQ_FATAL_TEST_MACRO( tokens.size() != m_no_columns,
			       0,
			       "uqDataSet::readFromFile",
			       "The number of data points differs from the number of columns in the header" );

	  for( unsigned int i = 0; i < tokens.size(); ++i )
	    {      
	      m_data_column[ m_column_name[i] ].push_back( StringToDouble( tokens[i] ) );
	    }

	  m_data_set_size++;

	}
    }
  
  data_in.close();

  m_file_name = fileName;

}

//********************************************************************
// uqDataSet::getValue
//********************************************************************
double uqDataSet::getValue( std::string nameVar, unsigned int idx )
{

  UQ_FATAL_TEST_MACRO( idx >= m_data_set_size,
		       0,
		       "uqDataSet::getValue",
		       "The index provided exceeds the size of the data set" );
  
  UQ_FATAL_TEST_MACRO( m_data_column.find( nameVar ) == m_data_column.end(),
		       0,
		       "uqDataSet::getValue",
		       "The variable inquired is not existent in the data set" );
  
  return m_data_column[ nameVar ][ idx ];

}

//********************************************************************
// uqDataSet::getMapDataPoint
//********************************************************************
void uqDataSet::getMapDataPoint( unsigned int idx, std::map<std::string,double>& map_data_point ) const
{

  map_data_point.clear();

  UQ_FATAL_TEST_MACRO( idx >= m_data_set_size,
		       0,
		       "uqDataSet::getMapDataPoint",
		       "The index provided exceeds the size of the data set" );

  for( std::map< std::string,std::vector<double> >::const_iterator it = m_data_column.begin();
       it != m_data_column.end(); ++it )
    {

      map_data_point[ it->first ] = it->second[ idx ];

    }

}


//********************************************************************
// uqDataSet::printDataSet
//********************************************************************
void uqDataSet::printInfo() const
{

  std::map< std::string,std::vector<double> >::const_iterator it; 

  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "File name: " << m_file_name << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "#id\t";
  for( unsigned int i = 0; i < m_no_columns; ++i )
    {
      std::cout << m_column_name[ i ] << "\t";
    }
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;

  for( unsigned int j = 0; j < m_data_set_size; ++j )
    {

      std::cout << (j+1) << "\t";
      for( unsigned int i = 0; i < m_no_columns; ++i )
	{
	  it = m_data_column.find( m_column_name[ i ] );
	  std::cout << it->second[ j ] << "\t";
	}
      std::cout << std::endl;

    }

  std::cout << "-------------------------------------------------------" << std::endl;

}

//********************************************************************
// uqDataSet::cleanContent
//********************************************************************
void uqDataSet::cleanContent()
{

  m_data_column.clear();
  m_column_name.clear();

  m_file_name = "";
  m_data_set_size = 0;
  m_no_columns = 0;

}

#endif // __UQ_DATASET_H
