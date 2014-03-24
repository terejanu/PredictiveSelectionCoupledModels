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
% generateData_oscillator.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear all;
close all;
clc;

%% ------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
randn('state',100);

%% ------------------------------------------------------------------------
% Measurement noise
%--------------------------------------------------------------------------
std_Q_oscillator = 0.1;

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
% Time when measurements are collected (obs are actual Forces)
%--------------------------------------------------------------------------
time_meas.t0 = pi/2; 
time_meas.dt = 1;
time_meas.tf = 4*pi - pi/2;
time_meas.tspan = time_meas.t0 : time_meas.dt : time_meas.tf;
time_meas.nMeas = length(time_meas.tspan);

%% ------------------------------------------------------------------------
% Hypotheses for the forcing term
%--------------------------------------------------------------------------
H_used = [1, 2, 3];

H{1}.name = 'modelOscillator_linearSpring';
H{1}.color = [1 0 0];
H{1}.no_param = 2;
H{1}.param_str = {'c', 'k10'};
H{1}.prior(1,1:2) = [0.05 0.15];                % c
H{1}.prior(2,1:2) = [3 5];                      % k10

H{2}.name = 'modelOscillator_cubicSpring';
H{2}.color = [0 1 0];
H{2}.no_param = 3;
H{2}.param_str = {'c', 'k30', 'k32'};
H{2}.prior(1,1:2) = [0.05 0.15];                % c
H{2}.prior(2,1:2) = [3 5];                      % k30
H{2}.prior(3,1:2) = [-6 -4];                    % k32

H{3}.name = 'modelOscillator_quinticSpring';
H{3}.color = [1 0 1];
H{3}.no_param = 4;
H{3}.param_str = {'c', 'k50', 'k52', 'k54'};
H{3}.prior(1,1:2) = [0.05 0.15];                % c
H{3}.prior(2,1:2) = [3 5 ];                     % k50
H{3}.prior(3,1:2) = [-6 -4];                    % k52
H{3}.prior(4,1:2) = [0.5 1.5];                  % k54

%% ------------------------------------------------------------------------
% Normalize priors
%--------------------------------------------------------------------------
for i = 1 : length(H)
    
    H{i}.normVal = diff(H{i}.prior,1,2);
    H{i}.priorNorm = H{i}.prior ./ repmat( H{i}.normVal, 1, 2 );
    
end;


%% ------------------------------------------------------------------------
% Generate the truth & measureements (OSCILLATOR)
%--------------------------------------------------------------------------
param_true = [0.1, 4, -5, 1];

IC = [0.5 0.5];
[t state_ev] = ode45(@(t,y) modelOscillator_quinticSpring(t,y,param_true),[0 time_meas.tspan],IC); % Solve ODE

f_meas = state_ev(:,2).^2/2;
% std_Q = f_meas * std_Q_oscillator;
% y_meas = f_meas + normrnd( 0, std_Q);
y_meas = f_meas .* lognrnd( 0, std_Q_oscillator, size(f_meas,1),size(f_meas,2) );

y_meas = y_meas(2:end);

% Plot Data (OSCILLATOR)

[t state_ev_quintic] = ode45(@(t,y) modelOscillator_quinticSpring(t,y,param_true),time.tspan,IC); % Solve ODE
f_quintic = state_ev_quintic(:,2).^2/2;

[t state_ev_linear] = ode45(@(t,y) modelOscillator_linearSpring(t,y,param_true(1:2)),time.tspan,IC); % Solve ODE
f_linear = state_ev_linear(:,2).^2/2;


figure; hold on;
plot( time.tspan, f_quintic, 'b' );
plot( time.tspan, f_linear, 'g' );
plot( time_meas.tspan, y_meas, '*r' );

figure; hold on;
plot( state_ev_quintic(:,1), state_ev_quintic(:,2), 'b' );
plot( state_ev_linear(:,1), state_ev_linear(:,2), 'g' );

%% ------------------------------------------------------------------------
% Save measurements into a file (OSCILLATOR)
%--------------------------------------------------------------------------
fid = fopen('oscillator_meas.txt', 'w');
fprintf(fid, 'scn_time obs_kinetic_energy\n');

for i = 1 : time_meas.nMeas
    fprintf(fid, '%f %f\n', time_meas.tspan(i), y_meas(i));
end;
    
fclose(fid);

