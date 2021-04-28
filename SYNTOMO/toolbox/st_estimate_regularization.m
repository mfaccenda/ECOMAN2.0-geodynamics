function [Ldamp,Lsmooth] = st_estimate_regularization(nobs,dt,long,latg,rzg,aPhase)
% ST_ESTIMATE_REGULARIZATION: Makes an estimate of the appropriate damping
% and smoothing multipliers for the regularization equations in SYNTOMO.
%
% INPUT
%     nobs: Number of travel-time observations
%       dt: Mean 1-sigma uncertainty of travel-time observations (s)
%     long: Longitudinal grid vector defining perturbational model (degrees)
%     latg: Latitudinal grid vector defining perturbational model (degrees)
%      rzg: Radial grid vector defining perturbational model (km)
%   aPhase: Define as 'P' or 'S' for compressional or shear wave speed
%           inversion.
%
% OUTPUT
%   Ldamp: Damping multiplier estimate
%   Lsmooth: Smoothing multiplier estimate
%
% NOTES
% + The estimate for Ldamp is the right order of magnitude but generally
%   over-estimated.
% + The estimate for Lsmooth is generally under-estimated by ~5.
%
% B. VanderBeek (OCT-2020)
%
% LICENSE
% Copyright (c) 2018-2020, Università di Padova, Manuele Faccenda
% All rights reserved.
%
% This software was developed at:
%   Dipartimento di Geoscienze
%   Università di Padova, Padova         
%   via Gradenigo 6,            
%   35131 Padova, Italy
%
% Project:    ECOMAN 
% Funded by:  ERC StG 758199 - NEWTON
% 
% ECOMAN is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, version 3 of the License.
% 
% ECOMAN is distributed WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details. You should have received a 
% copy of the GNU General Public License along with ECOMAN. If not, see 
% <http://www.gnu.org/licenses/>.
%
% Contact:
%   Manuele Faccenda    [manuele.faccenda@unipd.it]
%   Brandon VanderBeek  [brandonpaul.vanderbeek@unipd.it]
% 
% Main development team:
%   Manuele Faccenda    [manuele.faccenda@unipd.it]
%   Brandon VanderBeek  [brandon.vanderbeek@unipd.it]
%   Jianfeng Yang   
%   Albert de Montserrat Navarro
%

% Dimensions of perturbational array
n1     = length(long);
n2     = length(latg);
n3     = length(rzg);
nparam = n1*n2*n3;

% Estimate of average slowness perturbation required to fit data error
Re     = 6371;
kmpdeg = pi*Re/180;
dx     = kmpdeg*(max(long) - min(long))/(n1 - 1);
dy     = kmpdeg*(max(latg) - min(latg))/(n2 - 1);
dz     = (max(rzg) - min(rzg))/(n3 - 1);
dL     = sqrt((dx^2) + (dy^2) + (dz^2));
du     = dt/dL;
% Make a vector of random slowness perturbations with above standard deviation
du = du*randn(nparam,1);

% Get some reference velocities. Chosen 1D model does not significantly
% effect estimated regularization multipliers.
if strcmp(aPhase,'P')
    [vi,~] = get_TauP_model(abs(rzg),'prem');
elseif strcmp(aPhase,'S')
    [~,vi] = get_TauP_model(abs(rzg),'prem');
else
    error('Unrecognized value for ''aPhase''. Expected ''P'' or ''S''.');
end
vi = vi(:);
vi = permute(repmat(vi,1,n1,n2),[2,3,1]);
vi = vi(:);

% Damping Multiplier Estimate
% Solution damping prevents the inversion from fitting noise in the data
% and making perturbations in poorly sampled regions of the model. Because
% the inverse problem is weighted by the travel-time uncertainties, when
% the data is fit to the error the sum of the squared normalized residuals
% will equal the number of observations. To prevent over-fitting the data,
% a damping multiplier should be chosen such that random perturbations
% required to fit data noise result in an increase to the misfit function
% similar to the contribution from the data residuals.
Ldamp = sqrt(nobs/sum((du.*vi).^2)); % Damping multiplier. Note that damping is weighted by prior velocities

% Smoothing Multiplier Estimate
% The smoothing multipliers are chosen using the same logic as described
% for the damping multiplier.
% Construct Laplacian smoothing constraints
[Sx,isx,jsx] = LaplaceCoeff(long,1,[n1,n2,n3]);
[Sy,isy,jsy] = LaplaceCoeff(latg,2,[n1,n2,n3]);
[Sz,isz,jsz] = LaplaceCoeff(rzg,3,[n1,n2,n3]);
% Scale Laplacian coefficients
Sx = Sx.*vi(jsx);
Sy = Sy.*vi(jsy);
Sz = Sz.*vi(jsz);
% Build sparse smoothing matrices
Sx = sparse(isx,jsx,Sx,nparam,nparam);
Sy = sparse(isy,jsy,Sy,nparam,nparam);
Sz = sparse(isz,jsz,Sz,nparam,nparam);
% Define the perturbational vector roughness
Lsmooth = sum((cat(1,Sx,Sy,Sz)*du).^2);
% Estimated smoothing multiplier (constant in each dimension)
Lsmooth = sqrt(nobs/Lsmooth);

% Experience shows that the damping multiplier is over-estimated by a
% factor of ~5 while the smoothing operator is generally under-estimated by
% ~30. This is likely the result of using a random-valued vector as a 
% reference for the perturbational model.
Ldamp   = Ldamp/5;
Lsmooth = 30*Lsmooth;

% Only keep one significant figure in estimate
e       = floor(log10(Ldamp));
Ldamp   = round(Ldamp/(10^e))*(10^e);
e       = floor(log10(Lsmooth));
Lsmooth = round(Lsmooth/(10^e))*(10^e);
