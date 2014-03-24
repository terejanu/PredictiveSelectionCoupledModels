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
% test_convergence_chain_act_val.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear all;
close all;
clc;

addpath(genpath('mcmcdiag'));

my_chain = getChainFile( '../outputData/ForceSLD_rawChain_ml' );
run( '../outputData/ForceSLD_info' );

% Cummulative Plots
n_var = size(my_chain,2);

new_chain = zeros( size(my_chain) );

for i = 1 : n_var
    log_scale = params_min_max_log(i,3);
    norm_val =  params_min_max_log(i,2) -  params_min_max_log(i,1);
    if( log_scale == 1)
        new_chain(:,i) = 10.^( my_chain(:,i) * norm_val );
    else
        new_chain(:,i) = my_chain(:,i) * norm_val;
    end;
end;


for i = 1 : n_var
    plotCumPlotsChain(i, new_chain(:,i), 50);
end;
