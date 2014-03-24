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
% plotCumPlotsChain.m
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function plotCumPlotsChain(nvar, my_chain1, maxlag)

ac_chain1 = acorr(my_chain1,maxlag);

m1 = mean( my_chain1 );
m2 = mean( my_chain1.^2 );

cumsum1_m1 = cumsum( my_chain1 - m1 );
cumsum1_m2 = cumsum( my_chain1.^2 - m2 );


my_mu = mean(my_chain1);
my_std = std(my_chain1);

my_chain2 = my_mu + my_std * randn(size(my_chain1));

m1 = mean( my_chain2 );
m2 = mean( my_chain2.^2 );

cumsum2_m1 = cumsum( my_chain2 - m1 );
cumsum2_m2 = cumsum( my_chain2.^2 - m2 );

figure; 
    subplot(241); 
        [f,xi] = ksdensity(my_chain1);
        plot(xi,f/trapz(xi,f));
        title(['param #' num2str(nvar)]);
    subplot(242); 
        plot(my_chain1);
        title('my chain');
    subplot(243); 
        plot(cumsum1_m1);
        title('cumsum m1');
    subplot(244); 
        plot(cumsum1_m2);
        title('cumsum m2');
    subplot(245); 
        plot(ac_chain1); 
        title(['auto corr (lag=' num2str(maxlag) ')']);
    
    subplot(246); 
        plot(my_chain2);
        title('rnd chain');
    subplot(247); 
        plot(cumsum2_m1);
        title('cumsum m1 (rnd)');
    subplot(248); 
        plot(cumsum2_m2);
        title('cumsum m2 (rnd)');
    
    drawnow;