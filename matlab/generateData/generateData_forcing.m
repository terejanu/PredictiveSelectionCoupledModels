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
% generateData_forcing.m
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
std_Q_forcing = 0.1;

%% ------------------------------------------------------------------------
% Time for discretization (total time)
%--------------------------------------------------------------------------
time.t0 = 0; 
time.dt = .05;
time.tf = 30*pi;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% Time when measurements are collected (obs are actual Forces)
%--------------------------------------------------------------------------
% time_meas.t0 = pi; 
% time_meas.dt = pi/2;
% time_meas.tf = 3*pi;

time_meas.t0 = pi*3; 
time_meas.dt = pi/3;
time_meas.tf = 5*pi;

time_meas.tspan = time_meas.t0 : time_meas.dt : time_meas.tf;
time_meas.nMeas = length(time_meas.tspan);

%% ------------------------------------------------------------------------
% Hypotheses for the forcing term
%--------------------------------------------------------------------------
H_used = [1, 2, 3];

H{1}.name = 'modelForcing_simpleLinearDecay';
H{1}.color = [1 0 0];
H{1}.no_param = 2;
H{1}.param_str = {'F0', 'tau'};
H{1}.prior(1,1:2) = [0.5 2];                % F0
H{1}.prior(2,1:2) = [pi/3 4*pi];            % tau

H{2}.name = 'modelForcing_oscillatoryLinearDecay';
H{2}.color = [0 1 0];
H{2}.no_param = 4;
H{2}.param_str = {'F0', 'tau', 'alpha', 'omega'};
H{2}.prior(1,1:2) = [0.5 2];                % F0
H{2}.prior(2,1:2) = [pi/3 4*pi];            % tau
H{2}.prior(3,1:2) = [0.1 0.4];              % alpha
H{2}.prior(4,1:2) = [1 4];                  % omega

H{3}.name = 'modelForcing_simpleExponentialDecay';
H{3}.color = [0 0 1];
H{3}.no_param = 2;
H{3}.param_str = {'F0', 'tau'};
H{3}.prior(1,1:2) = [0.5 2];                % F0
H{3}.prior(2,1:2) = [pi/3 4*pi];            % tau

H{4}.name = 'modelForcing_oscillatoryExponentialDecay';
H{4}.color = [1 0 1];
H{4}.no_param = 4;
H{4}.param_str = {'F0', 'tau', 'alpha', 'omega'};
H{4}.prior(1,1:2) = [0.5 2];                % F0
H{4}.prior(2,1:2) = [pi/3 4*pi];            % tau
H{4}.prior(3,1:2) = [0.1 0.4];              % alpha
H{4}.prior(4,1:2) = [1 4];                  % omega

%% ------------------------------------------------------------------------
% Normalize priors
%--------------------------------------------------------------------------
for i = 1 : length(H)
    
    H{i}.normVal = diff(H{i}.prior,1,2);
    H{i}.priorNorm = H{i}.prior ./ repmat( H{i}.normVal, 1, 2 );
    
end;


%% ------------------------------------------------------------------------
% Generate the truth & measureements (FORCING)
%--------------------------------------------------------------------------
param_true = [1, 2*pi, 0.2, 2];

f_meas = modelForcing_oscillatoryExponentialDecay( time_meas.tspan, param_true );
exact_QoI = computeQoI_forcing( 'modelForcing_oscillatoryExponentialDecay', time.tspan, param_true );

y_meas = f_meas .* lognrnd( 0, std_Q_forcing, size(f_meas,1),size(f_meas,2) );


% Plot Data (FORCING)

f_all = modelForcing_oscillatoryExponentialDecay( time.tspan, param_true );

figure; hold on;
plot( time.tspan, f_all, 'b' );
plot( time_meas.tspan, y_meas, '*r' );


%% ------------------------------------------------------------------------
% Save measurements into a file (FORCING)
%--------------------------------------------------------------------------
fid = fopen('../../input_files/forcing_meas.txt', 'w');
fprintf(fid, 'scn_time obs_force\n');

for i = 1 : time_meas.nMeas
    fprintf(fid, '%f %f\n', time_meas.tspan(i), y_meas(i));
end;
    
fclose(fid);

