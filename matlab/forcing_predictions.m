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
% forcing_predictions.m
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
% Import chain  ---  SLD
%--------------------------------------------------------------------------
% my_chain = getChainFile( '../outputData/ForceSLD_rawChain_ml' );
% run( '../outputData/ForceSLD_info' );
% run( '../outputData/ForceSLD_fp_q_seq' );
% modelStr = 'modelForcing_simpleLinearDecay';
% saveStr = 'ForceSLD';

%% ------------------------------------------------------------------------
% Import chain  ---  SED
%--------------------------------------------------------------------------
% my_chain = getChainFile( '../outputData/ForceSED_rawChain_ml' );
% run( '../outputData/ForceSED_info' );
% run( '../outputData/ForceSED_fp_q_seq' );
% modelStr = 'modelForcing_simpleExponentialDecay';
% saveStr = 'ForceSED';

%% ------------------------------------------------------------------------
% Import chain  ---  OLD
%--------------------------------------------------------------------------
% my_chain = getChainFile( '../outputData/ForceOLD_rawChain_ml' );
% run( '../outputData/ForceOLD_info' );
% run( '../outputData/ForceOLD_fp_q_seq' );
% modelStr = 'modelForcing_oscillatoryLinearDecay';
% saveStr = 'ForceOLD';

%% ------------------------------------------------------------------------
% Import chain  ---  OED
%--------------------------------------------------------------------------
my_chain = getChainFile( '../outputData/ForceOED_rawChain_ml' );
run( '../outputData/ForceOED_info' );
run( '../outputData/ForceOED_fp_q_seq' );
modelStr = 'modelForcing_oscillatoryExponentialDecay';
saveStr = 'ForceOED';

%% ------------------------------------------------------------------------
% Time for discretization (total time)
%--------------------------------------------------------------------------
time.t0 = 0; 
time.dt = .01;
time.tf = 30*pi;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);

%% ------------------------------------------------------------------------
% Measurement noise
%--------------------------------------------------------------------------
std_Q_forcing = 0.1;
tresh_min_time_decay = 0.1;

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
param_true = [1, 2*pi, 0.2, 2];

f_truth = feval( 'modelForcing_oscillatoryExponentialDecay', time.tspan, param_true );
exact_QoI_intForce = computeQoI_forcing( 'modelForcing_oscillatoryExponentialDecay', time.tspan, param_true );

% compute min time decay
ind = find( f_truth <= tresh_min_time_decay );
exact_QoI_minTimeDecay = time.tspan( ind(1) );


%% ------------------------------------------------------------------------
% Import observations
%--------------------------------------------------------------------------

dataFile = importdata( '../input_files/forcing_meas.txt' );
my_obs = dataFile.data;


%% ------------------------------------------------------------------------
% Prediction observable
%--------------------------------------------------------------------------
qL = 0.005;
qU = 0.995;

qoiSmp_intForce = zeros(size(new_chain,1),1);
qoiSmp_minTimeDecay = zeros(size(new_chain,1),1);

funcSmp = zeros(size(new_chain,1),time.nSteps);
dataSmp = zeros(size(new_chain,1),time.nSteps);
qL_func = zeros(1, time.nSteps);
qU_func = zeros(1, time.nSteps);

qL_data = zeros(1, time.nSteps);
qU_data = zeros(1, time.nSteps);

for i = 1 : size(new_chain,1)
    param_smp = new_chain(i,:);
    std_smp = new_chain(i,end);
    
    funcSmp(i,:) = feval( modelStr, time.tspan, param_smp );
    dataSmp(i,:) = funcSmp(i,:) .* lognrnd( 0, sqrt(std_Q_forcing^2 + std_smp^2), ...
        size(funcSmp(i,:),1),size(funcSmp(i,:),2) );
    
    qoiSmp_intForce(i) = trapz( time.tspan, funcSmp(i,:) );
    
    % compute min time decay
    ind = find( funcSmp(i,:) <= tresh_min_time_decay );
    qoiSmp_minTimeDecay(i) = time.tspan( ind(1) );
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

    subplot(121); hold on;

        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,1));
        plot(xi,f,'r');
    
        [f,xi] = ksdensity(qoiSmp_intForce);
        plot(xi,f,'g');

        yl = ylim;
        line( [exact_QoI_intForce exact_QoI_intForce], yl );
        title('integrated force');
        
    subplot(122); hold on;
    
        [f,xi] = ksdensity(fp_mc_QoiSeq_unified(:,2));
        plot(xi,f,'r');

        [f,xi] = ksdensity(qoiSmp_minTimeDecay);
        plot(xi,f,'g');
        
        yl = ylim;
        line( [exact_QoI_minTimeDecay exact_QoI_minTimeDecay], yl );
        title('min time decay');        
        
saveas( gcf, ['figs/' saveStr '_qoi'], 'fig');