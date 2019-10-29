function PositionsX = chemotaxis_in_viscous_gradient(NUMBER_OF_CELLS,...
    SIMULATION_TIME, x, ViscosityGradient, NutrientsGradient)
% Simulates a population of Escherichia coli cells swimming in a 2D space
% with a gradient of viscosity and/or of alpha-methyl-aspartate concentration.
%
% 
%   Usage: 
%     PositionsX = chemotaxis_in_viscous_gradient(NUMBER_OF_CELLS,...
%       SIMULATION_TIME, x, ViscosityGradient, NutrientsGradient)
%  
%   Input:
%       NUMBER_OF_CELLS is the number of E. coli cells in the simulation.
%
%       SIMULATION_TIME is the duration in seconds of the simulation.
%
%       x is the spatial axis in micrometers accross the viscosity and/or
%           nutrients gradients.
%
%       ViscosityGradient is the viscosity in kg/(m s).along x.
%
%       NutrientsGradient is the alpha-methyl-aspartate concentration in mM.
%  
%   Output:
%     Positions_x is the final position of all bacteria at the end of the
%       simulation along x.
%
%   Calls: chemotactic_pathway_model , theoretical_friction_coefficients
%
%   Requires: Statistics toolbox.
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
% 

TimeStep = 0.1; % Time step in seconds of the simulation.
Iterations = round(SIMULATION_TIME/TimeStep);

%% Experimental constants for Escherichia coli AW405
CELL_LENGTH = 2; % Mean cell length in microns.
CELL_WIDTH = 0.75; % Mean cell width in microns.
TEMPERATURE = 33; % Temperature in celsius.
MEAN_RUN_LENGTH = 1; % Mean run length in seconds.
MEAN_TUMBLE_TIME = 0.2; % Mean tumble time in seconds.
MEAN_RUN_SPEED = 20; % Average speed during a run in water viscosity (micron/s).
MEAN_TUMBLE_ANGULAR_VELOCITY = 4.735; % Average angular velocity during a 
                                      % tumble in water (rad/s).

%% Gradients
% Reinterpolation of gradients to a 1µm grid
sx = min(x):max(x); % crossgradient dimension, in microns
sy = sx; % along gradient dimension, in microns
NutrientsGradient = interp1(x, NutrientsGradient, sx); % interpolates to equally spaced values
ViscosityGradient = interp1(x, ViscosityGradient, sx); % interpolates to equally spaced values

%% Behavioural gradients
% Gradients of swimming velocity and tumble angular velocity imposed by
% spatial changes in viscosity.

[ft, fr, ~, Dr]  = theoretical_friction_coefficients(...
    repmat(CELL_LENGTH, length(NutrientsGradient),1)*1e-6/2,...
    repmat(CELL_WIDTH, length(NutrientsGradient),1)*1e-6/2,...
    repmat(CELL_WIDTH, length(ViscosityGradient),1)*1e-6/2,...
    TEMPERATURE, 0, ViscosityGradient(:));

RunningSpeed = nan(size(ft,1),1);
TumbleAngularVelocity = nan(size(ft,1),1);

for ii=1:size(ft,1)
    RunningSpeed(ii) = sqrt(min(ft(:,2))./ft(ii,2)*MEAN_RUN_SPEED^2);
    TumbleAngularVelocity(ii) =...
        sqrt(min(fr(:,1))./fr(ii,1)*MEAN_TUMBLE_ANGULAR_VELOCITY^2);
end

%% Initial bacterial distributions
% Initial x and y positions and orientations of bacteria assuming uniform
% distributions.
Orientations = 2*pi*rand(NUMBER_OF_CELLS,1);
PositionsX = min(sx) + (range(sx))*rand(NUMBER_OF_CELLS,1);
Positions_y = min(sy) + (range(sy))*rand(NUMBER_OF_CELLS,1);

%% Methylation levels at steady state to start with
% Calculates methylation levels M at steady state for each location along
% across the nutrient gradient.
Methylation_level = ones(length(NutrientsGradient),1)*0.5;

for tt = 1:10000
    [Methylation_level,~,~] = chemotactic_pathway_model(Methylation_level,...
        NutrientsGradient, TimeStep,...
        MEAN_RUN_LENGTH, MEAN_TUMBLE_TIME);
end

%% Assigns the correct M level at each particle according to its position in x
[~,Locb] = ismember(NutrientsGradient(round(PositionsX)+1),NutrientsGradient);
Methylation_level = Methylation_level(Locb);

%% Change random number generator to 'simdTwister': SIMD-oriented Fast Mersenne Twister
rng(1,'simdTwister') % rng('default') %reverts change

%% Declare variables for IBM simulation
RunStarts = ones(NUMBER_OF_CELLS,1); % Vector of the iterations in which the
                                     % last run started for each cell.

%% IBM SIMULATION

for ii = 2:Iterations-1 % Main timestepping loop, starting at time = TimeStep.
    % Pitfall: if a cell executes 2 consecutive
    % tumbles, it will be recorded as 2 short tumbles
    % rather than as one long tumble.
    
    % Determines the methylation state of each cell and its probability of
    % tumbling.
    [Methylation_level,~,TumbleProbability] =...
        chemotactic_pathway_model(Methylation_level,...
        NutrientsGradient(ceil(PositionsX)), TimeStep,...
        MEAN_RUN_LENGTH, MEAN_TUMBLE_TIME);
    TumbleProbability(TumbleProbability>1) = 1;
    tumble = binornd(1,TumbleProbability); % Determines what cells start a 
                                           % tumble given their P (their
                                           % binomial probability of
                                           % tumbling).
    t = tumble(:) & RunStarts(:)<ii; % Index for the swimming cells that
                                     % started to tumble this time step
                                     % (and were not tumbling before).
    r = ~t(:) & RunStarts(:)<=ii; % Index for the cells that are running this
                                  % time step.
    
    %% Tumbles
    % Assign random tumble durations and orientations to cells that started
    % to tumble this time step according to the relevant distribution
    % (exponential for tumble duration and distribution and a normal.
    TumbleDurations = round(exprnd(MEAN_TUMBLE_TIME/TimeStep,sum(t),1));
    Orientations(t) = Orientations(t) +...
        (TumbleAngularVelocity(floor(PositionsX(t)+1)) .*TimeStep).*...
        randn(sum(t),1);
    
    %% Runs
    % Apply rotational brownian motion component and translational
    % displacement only to cells actively running.
    Orientations(r) = Orientations(r) +...
        sqrt(2*Dr(floor(PositionsX(r)+1),2) .* TimeStep) .* randn(sum(r),1);
    Positions_y(r) = Positions_y(r) +...
        RunningSpeed(floor(PositionsX(r)+1)) .* TimeStep .*sin(Orientations(r));
    PositionsX(r) = PositionsX(r) +...
        RunningSpeed(floor(PositionsX(r)+1)) .* TimeStep .*cos(Orientations(r));
    
    % Reverses running cells that bump into the x boundaries with the same
    % angle.
    ParticlesInXBoundaries = ...
        all([r, any([PositionsX < min(sx), (PositionsX  > max(sx))],2)], 2);
    Orientations(ParticlesInXBoundaries) = ...
        atan2(sin(Orientations(ParticlesInXBoundaries)),...
        - cos(Orientations(ParticlesInXBoundaries)));
    PositionsX(PositionsX<min(sx)) = ...
        min(sx) - PositionsX(PositionsX<min(sx));
    PositionsX(PositionsX>max(sx)) = ...
        max(sx)*2-PositionsX(PositionsX>max(sx));
    
    % Makes cells in upper or lower boundaries reappear in the opposite side.
    Positions_y(Positions_y>max(sy)) = ...
        Positions_y(Positions_y>max(sy))-range(sy);
    Positions_y(Positions_y<min(sy)) = ...
        Positions_y(Positions_y<min(sy))+range(sy);
    
    %% reset state of tumbling and bumped cells
    RunStarts(t) = ii + TumbleDurations + 1; % Reassign run_starts to cells that
                                             % have started tumbling.
end