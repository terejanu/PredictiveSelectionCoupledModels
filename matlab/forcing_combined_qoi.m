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
% forcing_combined_qoi.m
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
% TRUTH
%--------------------------------------------------------------------------
param_true = [1, 2*pi, 0.2, 2];

f_truth = feval( 'modelForcing_oscillatoryExponentialDecay', time.tspan, param_true );
exact_QoI_intForce = computeQoI_forcing( 'modelForcing_oscillatoryExponentialDecay', time.tspan, param_true );

% compute min time decay
ind = find( f_truth <= tresh_min_time_decay );
exact_QoI_minTimeDecay = time.tspan( ind(1) );

%% ------------------------------------------------------------------------
% Get individual qois
%--------------------------------------------------------------------------

% ------------ SLD
run( '../outputData/ForceSLD_fp_q_seq' );
SLD_qoi = fp_mc_QoiSeq_unified;

% ------------ SED
run( '../outputData/ForceSED_fp_q_seq' );
SED_qoi = fp_mc_QoiSeq_unified;

% ------------ OLD
run( '../outputData/ForceOLD_fp_q_seq' );
OLD_qoi = fp_mc_QoiSeq_unified;

% ------------ OED
% run( '../outputData/ForceOED_fp_q_seq' );
% OED_qoi = fp_mc_QoiSeq_unified;

%% ------------------------------------------------------------------------
% Get combined qois
%--------------------------------------------------------------------------

run( '../outputData/forcing_qoi' );
Combined_qoi = combined_qoi_seq;

%% ------------------------------------------------------------------------
% Plot pdfs
%--------------------------------------------------------------------------

figure; 
    
    subplot(121); hold on;

        [f,xi] = ksdensity(SLD_qoi(:,1));
        plot(xi,f,'b');

        [f,xi] = ksdensity(SED_qoi(:,1));
        plot(xi,f,'g');

        [f,xi] = ksdensity(OLD_qoi(:,1));
        plot(xi,f,'k');

%         [f,xi] = ksdensity(OED_qoi(:,1));
%         plot(xi,f,'m');

        [f,xi] = ksdensity(Combined_qoi(:,1));
        plot(xi,f,'r');

        yl = ylim;
        line( [exact_QoI_intForce exact_QoI_intForce], yl );
        title('integrated force');
        
        legend('SLD','SED','OLD','Combined');
%         legend('SLD','SED','OLD','OED','Combined');

    subplot(122); hold on;

        [f,xi] = ksdensity(SLD_qoi(:,2));
        plot(xi,f,'b');

        [f,xi] = ksdensity(SED_qoi(:,2));
        plot(xi,f,'g');

        [f,xi] = ksdensity(OLD_qoi(:,2));
        plot(xi,f,'k');

%         [f,xi] = ksdensity(OED_qoi(:,2));
%         plot(xi,f,'m');

        [f,xi] = ksdensity(Combined_qoi(:,2));
        plot(xi,f,'r');

        yl = ylim;
        line( [exact_QoI_minTimeDecay exact_QoI_minTimeDecay], yl );
        title('min time decay'); 
        
        legend('SLD','SED','OLD','Combined');
%         legend('SLD','SED','OLD','OED','Combined');