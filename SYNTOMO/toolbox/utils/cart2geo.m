function [LON,LAT,R] = cart2geo(X,Y,ELV,lon0,lat0)
% CART2GEO: Maps cartesian coordinates to geographic coordinates.
%
% INPUT
%      X: Array of cartesian x-coordinates to be mapped (km)
%      Y: Array of cartesian y-coordinates to be mapped (km)
%    ELV: Array of elevations (km)
%   lon0: Longitude of cartesian origin (deg.)
%   lat0: Latitude of cartesian origin (deg.)
%
% OUTPUT
%    LON: Longitude (deg.)
%    LAT: Latitude (deg.)
%      R: Radial coordinate (km)
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

% Earth Radius
Re = 6371;

% Get local z-coordinate
Z = sqrt(((Re+ELV).^2) - (X.^2) - (Y.^2));

% Store input array sizes
[ni,nj,nk] = size(Z);

% Convert cartesian coordinates to latitude/longitude
% If we place the Earth in a box, the X,Y-plane falls along the equator. In
% this global cartesian system, the local x, y, and z coordinates for a
% model centered on the equator/prime meridian correspond to ygloal,
% zglobal, and xglobal respectively. In this geometry latitude and
% longitude can be obtained from a cartesian to spherical coordinate
% transform. Assumes perfectly spherical Earth.

% Get global cartesian coordinates
xyz_global = rotz(lon0)*roty(-lat0)*[Z(:)';X(:)';Y(:)'];

% Convert to geographic coordinates
[LON,LAT,R] = cart2sph(xyz_global(1,:)',xyz_global(2,:)',xyz_global(3,:)');
LON = rad2deg(LON);
LAT = rad2deg(LAT);

% Return arrays in same dimension
LON = reshape(LON,ni,nj,nk);
LAT = reshape(LAT,ni,nj,nk);
R   = reshape(R,ni,nj,nk);
