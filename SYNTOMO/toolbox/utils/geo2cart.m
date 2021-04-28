function [X,Y,Z,ELEV] = geo2cart(LON,LAT,R,lon0,lat0,varargin)
% GEO2CART: Maps geographic coordinates to local cartesian coordinates.
%
% INPUT
%    LON: Array of longitude coordinates to be mapped (deg.)
%    LAT: Array of latitude coordinates to be mapped (deg.)
%      R: Radial coordinate
%   lon0: Longitude of cartesian origin (deg.)
%   lat0: Latitude of cartesian origin (deg.)
% <optional>
%      Re: Radius
%   alpha: Orientation of cartesian x-axis with respect to east
%          (CCW-positive; deg)
%
% OUTPUT
%      X: Cartesian x-coordinate
%      Y: Cartesian y-coordinate
%      Z: Cartesian z-coordinate
%   ELEV: Elevation with respect to surface of sphere
%
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

% This function returns the cartesian equivalent of the input geographic
% coordinates by placing the Earth in a box. The cartesian origin
% coresponds to (lon0,lat0,0).

if isempty(varargin)
    Re    = 6371;
    alpha = 0;
elseif length(varargin) == 1
    Re    = varargin{1};
    alpha = 0;
elseif length(varargin) == 2
    Re    = varargin{1};
    alpha = varargin{2};
else
    error('Incorrect number of variable input arguments.');
end

% Store input array sizes
[ni,nj,nk] = size(R);

% Convert to global cartesian coordinates
[Xglob,Yglob,Zglob] = sph2cart(deg2rad(LON),deg2rad(LAT),R);

% Get local cartesian coordinates
% Rotate origin to equator/prime meridian (Global cartesian)
xyz_local = (roty(-lat0)')*(rotz(lon0)')*[Xglob(:)';Yglob(:)';Zglob(:)']; 
% In global cartesian coordinates, x is into Earth, y is east, and z is
% north after the above rotation
xyz_local = xyz_local([2;3;1],:);

% Adjust for local coordinate system rotation
xyz_local = rotz(alpha)*xyz_local;

% Reshape back to input size
X = reshape(xyz_local(1,:),ni,nj,nk);
Y = reshape(xyz_local(2,:),ni,nj,nk);
Z = reshape(xyz_local(3,:),ni,nj,nk);

% Get elevation
ELEV = sqrt((X.^2) + (Y.^2) + (Z.^2)) - Re;

% Define Z as decreasing into box with surface at Z = 0
Z = Z - Re;
