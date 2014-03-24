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
% modelForcing_simpleLinearDecay.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function f = modelForcing_simpleLinearDecay(t,param)

%% ------------------------------------------------------------------------
% parameters
%--------------------------------------------------------------------------
F0 = param(1);
tau = param(2);

%% ------------------------------------------------------------------------
% forcing equation
%--------------------------------------------------------------------------
f = zeros(size(t));
if length(t) > 1,
    ind = find( t <= tau );
    f(ind) = F0*( 1-t(ind)/tau );
    ind = find( t > tau );
    f(ind) = 0;
else
    if t <= tf,
        f = F0*( 1-t/tau );
    else
        f = 0;
    end;
end;

