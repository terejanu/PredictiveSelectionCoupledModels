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
 * forcingCommon.h
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __FORCING_COMMON_H
#define __FORCING_COMMON_H

#define __FORCING_GSL_EPSABS                   1e-8
#define __FORCING_GSL_EPSREL                   1e-8
#define __FORCING_GSL_EPSROOT                  1e-3
#define __FORCING_QUADRATURE_INTERVALS         1000

#define __FORCING_INTEGRATION_PERIOD           100.0        // seconds
#define __FORCING_FORCE_DECAY_TOLERANCE        0.1         

#define __FORCING_MIN_SIGMA                    0.1          // known instrument error

#ifndef PI
#define PI 3.1415926535897932
#endif // PI

#endif // __FORCING_COMMON_H