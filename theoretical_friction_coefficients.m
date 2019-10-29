function [translationFrictionCoefficients,rotationalFrictionCoefficients,...
    translationalDiffusivities, rotationalDiffusivities] =...
    theoretical_friction_coefficients(...
    a, b, c, temperature, salinity, dinamycViscosity)

% Calculates translational and rotational friction coefficients, as well as
% translational and rotational diffusivities  using equations in appendix 6
% in Dusenbery 2009.
%
% Ex.:  [translationFrictionCoefficients,rotationalFrictionCoefficients,...
%           translationalDiffusivities, rotationalDiffusivities] =...
%           theoretical_friction_coefficients(...
%           a, b, c, temperature, salinity, dinamycViscosity)
% 
% Input variables:
%
%       a, b, c are the three semiaxes of the ellipsoid (in m)
%       temperature is temperature in Celsius degrees. Default = 33C.
%       salinity is salinity in 
%       dinamycViscosity is the dynamic viscosity in kg/m/s. 
%           Default = 7.4889e-04
%
% Output variables:
%
%       translationFrictionCoefficients is a (4X1) vector of the
%           translation friction coefficients parallel to each of the axes
%           and of the effective frictional coefficient calculated as the
%           harmonic mean of the coefficients of each of the three
%           orthogonal angles (Dusenbery 2009, Appendix 6).
%
%       rotationalFrictionCoefficients is a (3X1) vector of the rotational
%           friction coefficients about each of the axes.
%
%       translationalDiffusivities is a (4X1) vector of the translational
%           diffusivities plus the diffusion averaged over all orientations
%           (Dusenbery 2009, Appendix 6).
%
%       rotationalDiffusivities is a (3X1) vector of the rotational
%           diffusivities.
%
% Requires: SW_viscosity from Thermophysical properties of seawater toolbox
%        (http://web.mit.edu/seawater/)
%
% References
%   Dusenbery, D.B. (1998). Fitness landscapes for effects of
%       shape on chemotaxis and other behaviors of bacteria. Journal of
%       Bacteriology 180, 5978–5983.
%   Dusenbery, D.B. (2009). Living at microscale: the unexpected physics of
%       being small (Harvard University Press).
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
%  This file is part of trackbac, a collection of matlab scripts to
%  geometrically characterize and track swimming bacteria imaged by a
%  phase-contrast microscope

%% Constants
if nargin<4 || isempty(temperature)
    temperature = 33;
end
if nargin<5 || isempty(salinity)
    salinity = 0;
end
if nargin<6 || isempty(dinamycViscosity)
    dinamycViscosity = SW_Viscosity(temperature,'C',salinity,'ppt'); % kg/m-s
end

a = a(:);
b = b(:);
c = c(:);

BOLTZMANN_CONSTANT = 1.38e-23; % Boltzmann constant (J/K) or (m^2 kg s^-2 K^-1)
K = temperature+273.15; % temperature in Kelvins

%% Integrals for flow around ellipsoids (Dusenberry, appendix 6).
Sint = @(x) ((a.^2 + x).*(b.^2 + x).*(c.^2 + x)).^(-0.5);
salinity = integral(Sint,0,Inf,'ArrayValued',true);

r = [a, b, c];
Gint = @(x) ((repmat(a,1,3).^2 + x).*(repmat(b,1,3).^2 + x).*...
    (repmat(c,1,3).^2 + x)).^(-0.5) .* 1./(r.^2 + x);
G_E = r.^2 .* integral(Gint,0,Inf,'ArrayValued',true);

H_Ea = (G_E(:,2)+G_E(:,3))./(b.^2+c.^2); % Dusenberry 1998
H_Eb = (G_E(:,1)+G_E(:,3))./(a.^2+c.^2); % Dusenberry 1998
H_Ec = (G_E(:,1)+G_E(:,2))./(a.^2+b.^2); % Dusenberry 1998

%% Transational Friction Coefficients.
translationFrictionCoefficients(:,1) =...
    16*pi*dinamycViscosity./(salinity+G_E(:,1)); % Parallel to axis a eq A6.10.
translationFrictionCoefficients(:,2) =...
    16*pi*dinamycViscosity./(salinity+G_E(:,2)); % Parallel to axis b, eq A6.10.
translationFrictionCoefficients(:,3) =...
    16*pi*dinamycViscosity./(salinity+G_E(:,3)); % Parallel to axis c, eq A6.10.
translationFrictionCoefficients(:,4)...
    = 3./sum(1./translationFrictionCoefficients(:,1:3),2); % Harmonic mean.

%% Rotational friction coefficients.
rotationalFrictionCoefficients(:,1)=...
    16*pi*dinamycViscosity/3./H_Ea; % Friction constant for rotation about
                                    % the major semiaxis.
rotationalFrictionCoefficients(:,2)=...
    16*pi*dinamycViscosity/3./H_Eb; % Friction constant for rotation about the 
                                    % minor semiaxes.
rotationalFrictionCoefficients(:,3)=...
    16*pi*dinamycViscosity/3./H_Ec; % Friction constant for rotation about the
                                    % minor semiaxes.

%% Difusion coefficients
translationalDiffusivities =...
    BOLTZMANN_CONSTANT*K./translationFrictionCoefficients;
rotationalDiffusivities =...
    BOLTZMANN_CONSTANT*K./rotationalFrictionCoefficients;