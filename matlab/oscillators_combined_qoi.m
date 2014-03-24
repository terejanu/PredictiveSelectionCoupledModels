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
% oscillators_combined_qoi.m
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
% Time for discretization (total time)
%--------------------------------------------------------------------------
time.t0 = 0; 
time.dt = .05;
time.tf = 20*pi;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% Measurement noise
%--------------------------------------------------------------------------
std_Q_forcing = 0.1;
total_energy_eps = 0.1;

%% ------------------------------------------------------------------------
% TRUTH
%--------------------------------------------------------------------------
param_true = [0.1, 4, -5, 1];
IC = [0.5 0.5];

% [t state_ev] = ode45(@(t,y) modelOscillator_quinticSpring(t,y,param_true),[0 time_meas.tspan],IC); % Solve ODE
% f_truth = state_ev(:,2).^2/2;

% compute QoI
[t state_ev] = ode45(@(t,y) modelOscillator_quinticSpring(t,y,param_true), time.tspan,IC); % Solve ODE
kinetic = state_ev(:,2).^2/2;
total_energy = kinetic + potential_quinticSpring( state_ev(:,1), param_true );
ind = find( total_energy >= total_energy_eps );
exact_QoI_time_to_rest = time.tspan ( ind( end ) );

%% ------------------------------------------------------------------------
% Get individual qois
%--------------------------------------------------------------------------

% ------------ OLS
run( '../outputData/OscOLS_fp_q_seq' );
OLS_qoi = fp_mc_QoiSeq_unified;

% ------------ OCS
run( '../outputData/OscOCS_fp_q_seq' );
OCS_qoi = fp_mc_QoiSeq_unified;

% ------------ OQS
run( '../outputData/OscOQS_fp_q_seq' );
OQS_qoi = fp_mc_QoiSeq_unified;


%% ------------------------------------------------------------------------
% Get combined qois
%--------------------------------------------------------------------------

run( '../outputData/oscillator_qoi' );
Combined_qoi = combined_qoi_seq;

%% ------------------------------------------------------------------------
% Plot pdfs
%--------------------------------------------------------------------------

figure; hold on;

[f,xi] = ksdensity(OLS_qoi(:,1));
plot(xi,f,'b');

[f,xi] = ksdensity(OCS_qoi(:,1));
plot(xi,f,'g');

[f,xi] = ksdensity(OQS_qoi(:,1));
plot(xi,f,'k');

[f,xi] = ksdensity(Combined_qoi(:,1));
plot(xi,f,'r');

yl = ylim;
line( [exact_QoI_time_to_rest exact_QoI_time_to_rest], yl );
title('time to rest');
        
legend('OLS','OCS','OQS','Combined');
