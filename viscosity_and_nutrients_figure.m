% This is a simple script that loads nutrients and viscosity profiles and
% uses chemotaxis_in_viscous_gradient.m to simulate the behaviour of E.
% coli under different conditions.
%
%
%  Copyright (C) 2019,  Oscar Guadayol
%  oscar_at_guadayol.cat
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

load nutrients_and_viscosity % Loads the nutrients and viscosity profiles estimated using Comsol Multiphysics
NUMBER_OF_CELLS = 10000;
SIMULATION_TIME = 1200;
WaterViscosity = 0.001; % Water viscosity at 20C.
PosX_Nutrients_NoViscosity =...
    chemotaxis_in_viscous_gradient(NUMBER_OF_CELLS, SIMULATION_TIME,...
    distance, WaterViscosity*ones(size(distance)), nutrient_concentration1cp);
PosX_Nutrients_Viscosity =...
    chemotaxis_in_viscous_gradient(NUMBER_OF_CELLS, SIMULATION_TIME,...
    distance, viscosities, nutrient_concentration20cp);
PosX_NoNutrients_Viscosity =...
    chemotaxis_in_viscous_gradient(NUMBER_OF_CELLS, SIMULATION_TIME,...
    distance, viscosities, zeros(size(distance)));

NUMBER_OF_BINS = 100;
x_bins = min(distance):range(distance)/NUMBER_OF_BINS:max(distance);
MidXBins = x_bins +range(x_bins)/NUMBER_OF_BINS/2;
MidXBins = MidXBins(1:end-1);

%% Figure
fig = figure('Position', [736   568   764   423]);

ax(1) = axes; % Main plot with the bacterial abundance profiles.
c = get(ax(1),'colororder');
BackgroundAbundance = NUMBER_OF_CELLS/NUMBER_OF_BINS;
plot(ax(1), MidXBins, histcounts(PosX_Nutrients_NoViscosity, x_bins)./...
    BackgroundAbundance,'color',c(1,:))
hold(ax(1),'all')
plot(ax(1), MidXBins, histcounts(PosX_Nutrients_Viscosity,x_bins)./...
    BackgroundAbundance,'color',c(2,:))
plot(ax(1), MidXBins, histcounts(PosX_NoNutrients_Viscosity, x_bins)./...
    BackgroundAbundance,'color',c(3,:))
plot(ax(1), MidXBins, ones(size(MidXBins)),':k')

ax(1).XLim = [1 2000];
ax(1).XTickLabel = ax(1).XTick;
ax(1).YTickLabel = ax(1).YTick;
legend(ax(1),{'Nutrients no viscosity','Nutrients and viscosity',...
    'Viscosity gradient only'},'Location','NorthEast')
xlabel('Distance from source (\mum)')
ylabel([{'Bacterial abundance'}; {'(normalized to background)'}])
ax(1).FontSize = 11;
ax(1).Legend.FontSize = 11;
ax(1).Legend.Box = 'off';



ax(2) = axes('Position',[ 0.58 0.44 0.25 0.2]);  % Inset plot with the nutrients and viscosity profiles
ax(2).Box = 'on';
hold(ax(2),'on')
plot(ax(2), MidXBins, interp1(distance,nutrient_concentration20cp, MidXBins),'color',c(2,:))
plot(ax(2), MidXBins, interp1(distance,nutrient_concentration1cp, MidXBins),'color',c(1,:));
ax(3) = axes;
hold(ax(3),'on');
ax(3).Position = ax(2).Position;
plot(ax(3), MidXBins, interp1(distance,viscosities*1000, MidXBins),'color',c(4,:));
set(ax(3),'Color','none','YAxisLocation','right',...
    'YColor',c(4,:));
ax(2).XLim = [1 2000];
ax(3).XLim = ax(2).XLim;
ax(2).XTickLabel = ax(2).XTick;
ax(3).XTickLabel = '';
ax(2).YLabel.String = 'Nutrients (\muM)';
ax(3).YLabel.String = 'Viscosity (cP)';
ax(2).FontSize = 11;
ax(3).FontSize = 11;
