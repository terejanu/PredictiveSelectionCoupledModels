%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Copyright (C) 2011 Gabriel Terejanu - terejanu@cec.sc.edu
%
% This is an application framework for solving the problem of
% predictive model selection of coupled models as presented in the 
% following paper:
% 
% Gabriel Terejanu, Todd Oliver, Chris Simmons (2011). Application of 
% Predictive Model Selection to Coupled Models. In Proceedings of the World 
% Congress on Engineering and Computer Science 2011 Vol II, WCECS 2011, 
% pp. 927-932.
% 
% The framework is built on top of statistical library QUESO 
% (Quantification of Uncertainty for Estimation, Simulation and Optimization).
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
%--------------------------------------------------------------------------
%
% coupled_OQS_OED.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function dy = coupled_OQS_OED(t,y,param_OQS,param_OED)

global mass

%% ------------------------------------------------------------------------
% parameters
%--------------------------------------------------------------------------
c = param_OQS(1);
k50 = param_OQS(2);
k52 = param_OQS(3);
k54 = param_OQS(4);

m = mass;

%% ------------------------------------------------------------------------
% get force
%--------------------------------------------------------------------------
f = feval( 'modelForcing_oscillatoryExponentialDecay', t, param_OED );

%% ------------------------------------------------------------------------
% oscillator equation
%--------------------------------------------------------------------------
dy = zeros(2,1);
dy(1) = y(2);
dy(2) = f - c/m*y(2) - k50/m*y(1) - k52/m*y(1)^3 - k54/m*y(1)^5;

