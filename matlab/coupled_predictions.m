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
% coupled_predictions.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear all;
close all;
clc;

%% ------------------------------------------------------------------------
% Path
%--------------------------------------------------------------------------
addpath(genpath('mcmcdiag'));
addpath(genpath('generateData'));

%% ------------------------------------------------------------------------
% Object mass
%--------------------------------------------------------------------------
global mass
mass = 1;

%% ------------------------------------------------------------------------
% Import chain  ---  OLS - SED
%--------------------------------------------------------------------------
% my_chain_OLD = getChainFile( '../outputData/OscOLS_rawChain_ml' );
% my_chain_SED = getChainFile( '../outputData/ForceSED_rawChain_ml' );

run( '../outputData/OscOLS_ForceSED_fp_q_seq' );
saveStr = 'OscOLS_ForceSED';

%% ------------------------------------------------------------------------
% Import chain  ---  OLS - OLD
%--------------------------------------------------------------------------

% run( '../outputData/OscOLS_ForceOLD_fp_q_seq' );
% saveStr = 'OscOLS_ForceOLD';

%% ------------------------------------------------------------------------
% Import chain  ---  OCS - SED
%--------------------------------------------------------------------------

% run( '../outputData/OscOCS_ForceSED_fp_q_seq' );
% saveStr = 'OscOCS_ForceSED';

%% ------------------------------------------------------------------------
% Import chain  ---  OCS - OLD
%--------------------------------------------------------------------------

% run( '../outputData/OscOCS_ForceOLD_fp_q_seq' );
% saveStr = 'OscOCS_ForceOLD';

%% ------------------------------------------------------------------------
% Time for discretization (total time)
%--------------------------------------------------------------------------
time.t0 = 0; 
time.dt = .05;
time.tf = 30*pi;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% Measurement noise
%--------------------------------------------------------------------------
total_energy_eps = 0.1;

%% ------------------------------------------------------------------------
% TRUTH
%--------------------------------------------------------------------------
param_true_OQS = [0.1, 4, -5, 1];
param_true_OED = [1, 2*pi, 0.2, 2];
IC = [0.5 0.5];

% compute QoI
[t state_ev] = ode45(@(t,y) coupled_OQS_OED(t,y,param_true_OQS,param_true_OED), time.tspan, IC); % Solve ODE
kinetic = state_ev(:,2).^2/2;
total_energy = kinetic + potential_quinticSpring( state_ev(:,1), param_true_OQS );
ind = find( total_energy >= total_energy_eps );
exact_QoI_time_to_rest = time.tspan ( ind( end ) );
exact_QoI_max_displacement = max( abs( state_ev(:,1) ) );
exact_QoI_max_velocity = max( abs( state_ev(:,2) ) );

%% ------------------------------------------------------------------------
% Plot QoI
%--------------------------------------------------------------------------

figure; 

    subplot(131); hold on;    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,1));
        plot(xi,f,'r');
        yl = ylim;
        line( [exact_QoI_time_to_rest exact_QoI_time_to_rest], yl );
        title([saveStr ' : time to rest']);              
    subplot(132); hold on;    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,2));
        plot(xi,f,'r');
        yl = ylim;
        line( [exact_QoI_max_displacement exact_QoI_max_displacement], yl );
        title([saveStr ' : maximum displacement']);              
    subplot(133); hold on;    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,3));
        plot(xi,f,'r');
        yl = ylim;
        line( [exact_QoI_max_velocity exact_QoI_max_velocity], yl );
        title([saveStr ' : maximum velocity']);              
        
saveas( gcf, ['figs/' saveStr '_qoi'], 'fig');