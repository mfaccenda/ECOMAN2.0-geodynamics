function [vp,vss,vsf] = st_GetAnisoVel(model,xq,yq,zq,azim,elev)
% ST_GETANISOVEL: Computes seismic wave speeds at specified direction of
% propagation within an elastic model using the Christoffel equations.
% Created for SynTomo matlab package.
%
% INPUT
%   model: SynTomo model structure (see SynTomo manual)
%      xq: Vector of x-coordinates within model at which to compute wave 
%          speeds.
%      yq: Vector of y-coordinates within model at which to compute wave
%          speeds.
%      zq: Vector of z-coordinates within model at which to compute wave
%          speeds.
%    azim: Vector of propagation azimuths at which to compute wave speeds.
%          Azimuth is measured CCW from x-axis in DEGREES.
%    elev: Vector of propagation elevations at which to compute wave
%          speeds. Elevation is measured CCW from the x,y-plane in DEGREES.
%
% OUTPUT
%    vp: Compressional wave speed (km/s)
%   vss: Shear wave speed in slow polarization direction (km/s)
%   vsf: Shear wave speed in fast polarization direction (km/s)
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

% Interpolate elastic coefficients to query points
[~,w,ind] = TrilinInterp(model.long,model.latg,model.rzg,[],xq,yq,zq);
rho = sum(w.*model.RHO(ind),2);
c11 = sum(w.*model.C11(ind),2);
c12 = sum(w.*model.C12(ind),2);
c13 = sum(w.*model.C13(ind),2);
c14 = sum(w.*model.C14(ind),2);
c15 = sum(w.*model.C15(ind),2);
c16 = sum(w.*model.C16(ind),2);
c22 = sum(w.*model.C22(ind),2);
c23 = sum(w.*model.C23(ind),2);
c24 = sum(w.*model.C24(ind),2);
c25 = sum(w.*model.C25(ind),2);
c26 = sum(w.*model.C26(ind),2);
c33 = sum(w.*model.C33(ind),2);
c34 = sum(w.*model.C34(ind),2);
c35 = sum(w.*model.C35(ind),2);
c36 = sum(w.*model.C36(ind),2);
c44 = sum(w.*model.C44(ind),2);
c45 = sum(w.*model.C45(ind),2);
c46 = sum(w.*model.C46(ind),2);
c55 = sum(w.*model.C55(ind),2);
c56 = sum(w.*model.C56(ind),2);
c66 = sum(w.*model.C66(ind),2);

% Compute apparent velocities for each query point
npt = length(xq);
vp  = zeros(npt,1);
vss = zeros(npt,1);
vsf = zeros(npt,1);
for ii = 1:npt
    % Elastic tensor at query point
    Cij = [c11(ii),c12(ii),c13(ii),c14(ii),c15(ii),c16(ii);...
        c12(ii),c22(ii),c23(ii),c24(ii),c25(ii),c26(ii);...
        c13(ii),c23(ii),c33(ii),c34(ii),c35(ii),c36(ii);...
        c14(ii),c24(ii),c34(ii),c44(ii),c45(ii),c46(ii);...
        c15(ii),c25(ii),c35(ii),c45(ii),c55(ii),c56(ii);...
        c16(ii),c26(ii),c36(ii),c46(ii),c56(ii),c66(ii)];
    
    % Solve Christoffel equations
    [a,bs,bf] = christoffel(Cij,rho(ii),[azim(ii),elev(ii)].*(180/pi),false);
    vp(ii)    = a;
    vss(ii)   = bs;
    vsf(ii)   = bf;
end
