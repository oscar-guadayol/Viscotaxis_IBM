function [methylationLevel, activity, tumblingProbability] =...
    chemotactic_pathway_model(methylationLevel, aspartateConcentration,...
    timeStep, meanRunLength, meanTumbleTime)

% Model of the Escherichia coli chemosensory circuit.
%   Usage: 
%     [methylationLevel, activity, tumblingProbability] =...
%       chemotactic_pathway_model(methylationLevel, aspartateConcentration,...
%       TIME_STEP, MEAN_RUN_LENGTH, MEAN_TUMBLE_TIME)
%
% Input variables:
%
%     methylationLevel is the average methylation level of the receptor.
%     activity is average activity of the receptor's kinase complex.
%     aspartateConcentration is the alpha-methyl-aspartate concentration in
%         microM.
%     timeStep is the time step in seconds.
%     meanRunLength is the mean run length in seconds at the adapted state
%         (i.e. homogeneous nutrient field).
%     meanTumbleTime is the mean tumble time in seconds at the adapted state
%         (i.e. homogeneous nutrient field).
%  
% Output variables:
%
%     methylationLevel is the average methylation level of the receptor.
%     activity is average activity of the receptor's kinase complex.
%     tumblingProbability is the probability of tumbling during this time step.
%
%  References
%
% Kalinin, Y.V., Jiang, L., Tu, Y., and Wu, M. (2009). Logarithmic Sensing
%   in Escherichia coli Bacterial Chemotaxis. 96, 2439–2448.
%
% Shimizu, T.S., Tu, Y., and Berg, H.C. (2010). A modular gradient‐sensing
%   network for chemotaxis in Escherichia coli revealed by responses to
%   time‐varying stimuli. Molecular Systems Biology 6, 382.
%
% Tu, Y., Shimizu, T.S., and Berg, H.C. (2008). Modeling the chemotactic
%   response of Escherichia coli to time-varying stimuli. PNAS 105, 14855–14860.
%
% Vladimirov, N., Løvdok, L., Lebiedz, D., and Sourjik, V. (2008).
%   Dependence of Bacterial Chemotaxis on Gradient Shape and Adaptation
%   Rate. PLOS Computational Biology 4, e1000242.
%
% Cluzel, P., Surette, M., and Leibler, S. (2000). An Ultrasensitive
%   Bacterial Motor Revealed by Monitoring Signaling Proteins in Single
%   Cells. Science 287, 1652–1655.
%
%  Copyright (C) 2019,  Oscar Guadayol
%  oscar_at_guadayol.cat
%
%
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License, version 3.0, as
%  published by the Free Software Foundation.
% 
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License along
%  with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Chemotactic pathway constants
% From Kalinin et al 2009, Tu et al 2008, Shimizu et al 2010.

% Phosporilation constants %
N_RECEPTORS = 6; % number of responding receptor dimers
                 % in the MCP (Methyl-accepting 
                 % Chemotactic Proteins)
K_I = 18; % Dissociation constant of the ligand to the
          % inactive receptor in microM.
K_A = K_I/0.0062; % K_A = K_I/C dissociation constant of the
                  % ligand to the active receptor
alpha = 1.7; % Is the free-energy change per added methyl group (in units of kT).
             % In Shimizu et al alpha = 2kT at ~22C.
m_0 = 1; % methylation level at which the methylation level dependent 
         % freeEnergyDifference is 0.

% Methylation constants %
k_R = 0.005; % s^-1; Methylation rate for the inactive receptors.
k_B = 0.005; % s^-1; De-methylation rate for the active receptors.

% Tumble probability constants %
H = 10.3; % hill coefficient;
a_0 = 1/2; % The steady-state activity in the absence of stimuli.
           % In Shimizu et al a_0 = 1/3 at 32C.
           
%% Phosphorilation
freeEnergyDifference = alpha*(m_0 - methylationLevel); % Methylation level
                                                       % dependent free
                                                       % energy difference.
freeEnergy = N_RECEPTORS*(freeEnergyDifference +...
    log(1+aspartateConcentration(:)/K_I)-log(1+aspartateConcentration(:)/K_A));
activity = (1+exp(freeEnergy)).^-1; % Average kinase activity.

%% Methylation
% Linear model (Kalinin et al 2009).
F = k_R*(1-activity)-k_B*activity; % Function of dependence of methylatiol level
                                   % to receptor-kinase output.

% Michaelis-Menten model (Tu et al 2013, Shimizu et al 2010) at 32C.
% r_B = 2.7;
% a_B = 0.74;
% K_R = 0.43;
% % K_R = 0.2; % Value to make the function fit the one plotted in Shimizu et al 2010 fig 5B
% K_B = 0.30;
% V_R = 0.030;
% % V_R = 0.0357;% Value to make the function fit the one plotted in Shimizu et al 2010 fig 5B
% V_B0 = 0.030;
% phi = zeros(size(a));
% phi(a-a_B>0) = 1; % unit step function
% V_B = V_B0*(1+phi.*(a-a_B)./(1-a_B)*r_B);
% F = V_R*(1-a)./(K_R+1-a) - V_B.*a./(K_B+a); % The output of this function does not exactly match the function plotted in Shimizu et al 201 Fig 5B

methylationLevel = methylationLevel + F*timeStep; % Average methylation level.

%% Tumble probability (Cluzel et al 2000)

cwBiasHomogeneous = meanTumbleTime/(meanTumbleTime+meanRunLength); 
                  % Clockwise bias in the absence of a chemoattactant gradient.

aHalf = a_0*(1/cwBiasHomogeneous-1)^(1/H); % Kinase activity that induces a
                                           % clockwise bias of 0.5 in the
                                           % absence of a chemoattactant
                                           % gradient.
                                           
cwBias = 1./((aHalf./activity).^H+1); % Probability state for a cell of going
                                      % into tumbling.
                                      
runLength = meanTumbleTime./cwBias-meanTumbleTime; % Expected run length
                                                       % for each cell at each
                                                       % time step in seconds.
                                                       
frequency = 1./runLength; % Frequency of the Poisson process of run-tumble
                          % (Vladimirov et al 2009).
                          
tumblingProbability = timeStep.*frequency; % Probability of tumbling during
                                            % this time step.
