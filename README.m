% This is a set of scripts and function that model the behaviour of a
% population of Escherichia coli cells in the presence of a gradient of
% viscosity and/or of nutrients concentration. The main function is
% chemotaxis_in_viscous_gradient, which simulates population of individual
% E. coli cells swimming in a 2D space with a gradient of viscosity and/or
% of alpha-methyl-aspartate concentration. The response of individual cells
% to the nutrients is modelled in chemotactic_pathway_model, that codifies
% the E. coli's chemosensory circuit. the effect of varying viscosity is
% simulated with the ellipsoidal model for rotational and translational
% friction coefficients (theoretical_friction_coefficients).
% an example usage of these scripts is presented in 
%
%  Copyright (C) 2019,  Oscar Guadayol <oscar_at_guadayol.cat>
%
% LICENSE:
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


