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
% oscillators_predictions.m
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
% Import chain  ---  OLS
%--------------------------------------------------------------------------
my_chain = getChainFile( '../outputData/OscOLS_rawChain_ml' );
run( '../outputData/OscOLS_info' );
run( '../outputData/OscOLS_fp_q_seq' );
modelStr = 'modelOscillator_linearSpring';
saveStr = 'OscOLS';

%% ------------------------------------------------------------------------
% Import chain  ---  OCS
%--------------------------------------------------------------------------
% my_chain = getChainFile( '../outputData/OscOCS_rawChain_ml' );
% run( '../outputData/OscOCS_info' );
% run( '../outputData/OscOCS_fp_q_seq' );
% modelStr = 'modelOscillator_cubicSpring';
% saveStr = 'OscOCS';

%% ------------------------------------------------------------------------
% Import chain  ---  OQS
%--------------------------------------------------------------------------
% my_chain = getChainFile( '../outputData/OscOQS_rawChain_ml' );
% run( '../outputData/OscOQS_info' );
% run( '../outputData/OscOQS_fp_q_seq' );
% modelStr = 'modelOscillator_quinticSpring';
% saveStr = 'OscOQS';

%% ------------------------------------------------------------------------
% Time for discretization (total time)
%--------------------------------------------------------------------------
time.t0 = 0; 
time.dt = .05;
time.tf = 100;
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
% Measurement noise
%--------------------------------------------------------------------------
std_Q_oscillator = 0.1;
total_energy_eps = 0.1;

%% ------------------------------------------------------------------------
% Cummulative Plots
%--------------------------------------------------------------------------
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
    saveas( gcf, ['figs/' saveStr '_param_' num2str(i)], 'fig');
end;

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
exact_QoI_max_displacement = max( abs( state_ev(:,1) ) );
exact_QoI_max_velocity = max( abs( state_ev(:,2) ) );

%% ------------------------------------------------------------------------
% Import observations
%--------------------------------------------------------------------------
dataFile = importdata( '../input_files/oscillator_meas.txt' );
my_obs = dataFile.data;

%% ------------------------------------------------------------------------
% Prediction observable
%--------------------------------------------------------------------------
qL = 0.005;
qU = 0.995;

qoiSmp_time_to_rest = zeros(size(new_chain,1),1);
qoiSmp_max_displacement = zeros(size(new_chain,1),1);
qoiSmp_max_velocity = zeros(size(new_chain,1),1);
qoiSmp_max_velocity_model_unc = zeros(size(new_chain,1),1);

funcSmp = zeros(size(new_chain,1),time.nSteps);
dataSmp = zeros(size(new_chain,1),time.nSteps);
qL_func = zeros(1, time.nSteps);
qU_func = zeros(1, time.nSteps);

qL_data = zeros(1, time.nSteps);
qU_data = zeros(1, time.nSteps);

for i = 1 : size(new_chain,1)
    i
    param_smp = new_chain(i,:);
    std_smp = new_chain(i,end);
    
    % change this
    switch saveStr
        case 'OscOLS'
            [t state_ev] = ode45(@(t,y) modelOscillator_linearSpring(t,y,param_smp), time.tspan,IC); % Solve ODE
        case 'OscOCS'
            [t state_ev] = ode45(@(t,y) modelOscillator_cubicSpring(t,y,param_smp), time.tspan,IC); % Solve ODE
        case 'OscOQS'
            [t state_ev] = ode45(@(t,y) modelOscillator_quinticSpring(t,y,param_smp), time.tspan,IC); % Solve ODE
    end;
    kinetic_energy = state_ev(:,2).^2/2;
    funcSmp(i,:) = kinetic_energy;
    
    dataSmp(i,:) = funcSmp(i,:) .* lognrnd( 0, sqrt(std_Q_oscillator^2 + std_smp^2), ...
        size(funcSmp(i,:),1),size(funcSmp(i,:),2) );
    
    dataSmp_model_unc = funcSmp(i,:) .* lognrnd( 0, std_smp, ...
        size(funcSmp(i,:),1),size(funcSmp(i,:),2) );
    
    % compute QoI
    switch saveStr
        case 'OscOLS'
            potential_energy = potential_linearSpring( state_ev(:,1), param_smp );
        case 'OscOCS'
            potential_energy = potential_cubicSpring( state_ev(:,1), param_smp );
        case 'OscOQS'
            potential_energy = potential_quinticSpring( state_ev(:,1), param_smp );
    end;
    total_energy = kinetic_energy + potential_energy;
    ind = find( total_energy >= total_energy_eps );
    qoiSmp_time_to_rest(i) = time.tspan ( ind( end ) );
    qoiSmp_max_displacement(i) = max( abs( state_ev(:,1) ) );
    qoiSmp_max_velocity(i) = max( abs( state_ev(:,2) ) );
    qoiSmp_max_velocity_model_unc(i) = max( abs( sqrt(2*dataSmp_model_unc) ) );
end;

% get mean function
meanFunc = mean(funcSmp);

% get quantile functions
for i = 1 : time.nSteps
    qL_func(i) = quantile( funcSmp(:,i), qL );
    qU_func(i) = quantile( funcSmp(:,i), qU );
    qL_data(i) = quantile( dataSmp(:,i), qL );
    qU_data(i) = quantile( dataSmp(:,i), qU );
end;

% plot observable
figure; hold on;
plot( time.tspan, meanFunc, 'g');
plot( time.tspan, qL_func, '-.r', time.tspan, qU_func, '-.r');
plot( time.tspan, qL_data, '-.b', time.tspan, qU_data, '-.b');
plot( my_obs(:,1), my_obs(:,2), '*' );
hold off;
title('observable');
saveas( gcf, ['figs/' saveStr '_obs'], 'fig');

% plot qois
figure; 

    subplot(131); hold on;    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,1));
        plot(xi,f,'r');
        [f,xi] = ksdensity(qoiSmp_time_to_rest);
        plot(xi,f,'g');
        yl = ylim;
        line( [exact_QoI_time_to_rest exact_QoI_time_to_rest], yl );
        title([saveStr ' : time to rest']);              
    subplot(132); hold on;    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,2));
        plot(xi,f,'r');
        [f,xi] = ksdensity(qoiSmp_max_displacement);
        plot(xi,f,'g');
        yl = ylim;
        line( [exact_QoI_max_displacement exact_QoI_max_displacement], yl );
        title([saveStr ' : maximum displacement']);              
    subplot(133); hold on;    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,3));
        plot(xi,f,'r');
        [f,xi] = ksdensity(qoiSmp_max_velocity);
        plot(xi,f,'g');
        [f,xi] = ksdensity(qoiSmp_max_velocity_model_unc);
        plot(xi,f,'m');        
        yl = ylim;
        line( [exact_QoI_max_velocity exact_QoI_max_velocity], yl );
        title([saveStr ' : maximum velocity']);         

saveas( gcf, ['figs/' saveStr '_qoi'], 'fig');