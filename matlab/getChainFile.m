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
% getChainFile.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function myChain = getChainFile( fileName )

% load the data in the memory
run( fileName );

% get the chain from the last level
a = who('-regexp','rawChain_unified');

all_levels = [];
s1 = regexp(a, '_');

for i = 1 : length(s1)
    all_levels = [all_levels str2num(a{i}( s1{i}(2)+1 : s1{i}(3)-1 )) ];
end;

a = who('-regexp',['ip_ml_' num2str(max(all_levels)) '_rawChain_unified']);
myChain = eval(a{1});