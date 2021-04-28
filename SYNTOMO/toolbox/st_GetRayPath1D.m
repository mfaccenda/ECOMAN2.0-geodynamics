function ray = st_GetRayPath1D(aModel,aPhase,slat,slon,sdpt,rlat,rlon,latlims,lonlims,zlims,dsamp,Re)
% ST_GETRAYPATH1D: Defines ray path coordinates within local model. Ray path
% is computed via TauP assuming a 1D radial Earth model. Created for
% SynTomo Matlab package.
%
% INPUT
%    aModel: Name of 1D radial Earth model (see TauP manual for options)
%    aPhase: A seismic phase name (see TauP manual for options)
%      slat: Source latitude (deg.)
%      slon: Source longitude (deg.)
%      sdpt: Source depth (km)
%      rlat: Receiver latitude (deg.)
%      rlon: Receiver longitude (deg.)
%   latlims: A 2-component vector specifying minimum and maximum latitude
%            of model grid (deg.)
%   lonlims: A 2-component vector specifying minimum and maximum longitude
%            of model grid (deg.)
%     zlims: A 2-component vector specifying minimum and maximum
%            z-coordinate of model grid (deg.)
%     dsamp: A constant along-ray sampling interval (km)
%        Re: Earth radius (km)
%
% OUTPUT
%   ray: Structure containing ray path information with the following
%        fields...
%        tt_taup: The TauP predicted travel-time through 'aModel'
%              L: Length of ray path from last surface point (km)
%            lon: Longitude of ray (deg.)
%            lat: Latitude of ray (deg.)
%             rz: Radial ray coordinate (km)
%              x: Cartesian x-coordinate of ray (km)
%              y: Cartesian y-coordinate of ray (km)
%              z: Cartesian z-coordinate of ray (km)
%              d: Distance along ray from source (km)
%            azm: Ray azimuth (rad.)
%            elv: Ray elevation (rad.)
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

% Get 1D ray path via TauP
[rlat,rlon,rzray,ttray] = taup_path(aModel,aPhase,slat,slon,sdpt,rlat,rlon);
% Store travel-time
ray.tt_taup = ttray(end);

% Ray path in global cartesian coordinates
[xray,yray,zray] = geo2cart(rlon,rlat,Re + rzray,0,0,Re,0);
% Total ray length (needed for finite frequency kernels). We count from
% last surface point.
i0    = find(rzray(1:end-1) == 0,1,'last');
if isempty(i0)
    i0 = 1;
end
ray.L = sum(sqrt((diff(xray(i0:end)).^2) + (diff(yray(i0:end)).^2) + (diff(zray(i0:end)).^2)));

% Select all ray nodes within model
isin = ((rlon >= lonlims(1)) & (rlon <= lonlims(2)))...
    & ((rlat >= latlims(1)) & (rlat <= latlims(2)))...
    & ((rzray >= zlims(1)) & (rzray <= zlims(2)));

% Find the intersection of the ray path with the model boundary
iend = find(isin,1,'first');
nray = [(rlon(iend - 1) - rlon(iend)); (rlat(iend - 1) - rlat(iend));...
    rzray(iend - 1) - rzray(iend)];
[loni,lati,rzi] = get_box_intersection(nray,[rlon(iend),rlat(iend),rzray(iend)],...
    lonlims(1),lonlims(2),latlims(1),latlims(2),zlims(1),zlims(2));
[xi,yi,zi] = geo2cart(loni,lati,Re + rzi,0,0,Re,0);

% Subset ray path
ray.lon  = cat(1,loni,rlon(isin));
ray.lat  = cat(1,lati,rlat(isin));
ray.rz   = cat(1,rzi,rzray(isin));
ray.x    = cat(1,xi,xray(isin));
ray.y    = cat(1,yi,yray(isin));
ray.z    = cat(1,zi,zray(isin));
% Distance along ray path
ray.d    = sqrt((diff(ray.x).^2) + (diff(ray.y).^2) + (diff(ray.z).^2));
ray.d    = cumsum(ray.d);
ray.d    = cat(1,0,ray.d);

% Re-sample ray path
% Interpolate to a constant along-ray sampling interval
nq      = 1 + floor((ray.d(end) - ray.d(1))/dsamp);
dq      = linspace(ray.d(1),ray.d(end),nq)';
ray.lon = interp1(ray.d,ray.lon,dq,'linear');
ray.lat = interp1(ray.d,ray.lat,dq,'linear');
ray.rz  = interp1(ray.d,ray.rz,dq,'linear');
ray.x   = interp1(ray.d,ray.x,dq,'linear');
ray.y   = interp1(ray.d,ray.y,dq,'linear');
ray.z   = interp1(ray.d,ray.z,dq,'linear');
ray.d   = dq;
% Define ray orientations at segment nodes
dx = [ray.x(2)-ray.x(1); ray.x(3:end) - ray.x(1:end-2); ray.x(end)-ray.x(end-1)];
dy = [ray.y(2)-ray.y(1); ray.y(3:end) - ray.y(1:end-2); ray.y(end)-ray.y(end-1)];
dz = [ray.z(2)-ray.z(1); ray.z(3:end) - ray.z(1:end-2); ray.z(end)-ray.z(end-1)];
dL = sqrt((dx.^2) + (dy.^2) + (dz.^2));
dx = dx./dL;
dy = dy./dL;
dz = dz./dL;
[ray.azm,ray.elv] = cart2sph(dx,dy,dz);
