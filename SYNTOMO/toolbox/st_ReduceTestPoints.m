function [XG,YG,ZG,IND] = st_ReduceTestPoints(T,ray,vray,minx,maxx,miny,maxy,Re,XG,YG,ZG)
% ST_REDUCETESTPOINTS: Selects a subset of model points to consider when
% defining the HFFK in order to reduce computational time. Subset is 
% defined by a cone around ray path. Created for SynTomo Matlab package.
%
% INPUT
%      T: Wave period used to defin cone radius (s)
%    ray: SynTomo raypath structure (see function 'st_GetRayPath1D')
%   vray: Wave speeds along ray path (km/s)
%   minx: Minimum x-coordinate of model array
%   maxx: Maximum x-coordinate of model array
%   miny: Minimum y-coordinate of array
%   maxy: Maximum y-coordinate of array
%     Re: Earth radius
%     XG: Cartesian x-coordinate mesh
%     YG: Cartesian y-coordinate mesh
%     ZG: Cartesian z-coordinate mesh
%
% OUTPUT
%    XG: Subset of model x grid points
%    YG: Subset of model y grid points
%    ZG: Subset of model z grid points
%   IND: Linear index of node subset in mesh
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

% Define grid parameters. Note that cartesian grid-spacing actually varies
% with position because the model is defined in spherical coordinates but
% we approximate it as constant for selecting nodes.
[nx,ny,nz] = size(XG);
dlon       = (maxx - minx)./(nx - 1);
dlat       = (maxy - miny)./(ny - 1);
dx         = Re*pi*dlon/180;
dy         = Re*pi*dlat/180;

% Define maximum and minimum Fresnel radii
lray   = ray.d(end) - ray.d;
Rfx    = sqrt(T*vray(1).*lray(1).*(ray.L - lray(1))./ray.L)/sqrt(2);
Rfn    = sqrt(T*vray(end-1).*lray(end-1).*(ray.L - lray(end-1))./ray.L)/sqrt(2);
Rfx    = max(Rfx,(dx + dy)/2);
Rfn    = max(Rfn,(dx + dy)/2);

% First subset: simple box around ray path that encompasses maximum and
% minimum Fresnel radii.
% Nodal width of minimum Fresnel radius 
ni1 = ceil(Rfn/dx);
nj1 = ceil(Rfn/dy);
% Nodal width of maximum Fresnel radius 
ni2 = ceil(Rfx/dx);
nj2 = ceil(Rfx/dy);
% Index nearest station
ista = 1 + round((ray.lon(end) - minx)/dlon);
jsta = 1 + round((ray.lat(end) - miny)/dlat);
% Index nearest ray entry point
iin  = 1 + round((ray.lon(1) - minx)/dlon);
jin  = 1 + round((ray.lat(1) - miny)/dlat);
% Indices defining box
ii = [min(max(ista-ni1,1),max(iin-ni2,1)),max(min(ista+ni1,nx),min(iin+ni2,nx))];
jj = [min(max(jsta-nj1,1),max(jin-nj2,1)),max(min(jsta+nj1,ny),min(jin+nj2,ny))];
% Subset points
XG = XG(ii(1):ii(2),jj(1):jj(2),:); XG = XG(:);
YG = YG(ii(1):ii(2),jj(1):jj(2),:); YG = YG(:);
ZG = ZG(ii(1):ii(2),jj(1):jj(2),:); ZG = ZG(:);
% Define nodel index in model of subset
[JY,IX,KZ] = meshgrid(jj(1):jj(2),ii(1):ii(2),1:nz);
IND        = sub2ind([nx,ny,nz],IX,JY,KZ);

% Second subset: conical region around ray path
% Cone vertex vector (the line connecting the station and entry point)
n = -[ray.x(end)-ray.x(1); ray.y(end)-ray.y(1); ray.z(end)-ray.z(1)];
% Cone angle
theta = atan(Rfx/norm(n));
% Make unit vector
n = n./norm(n);
% Shift cone origin such that it coincides with the station and that the
% cone radius at the station is equal to the minimum Fresnel radius.
dL = Rfn/tan(theta);
x0 = ray.x(end) - dL*n(1);
y0 = ray.y(end) - dL*n(2);
z0 = ray.z(end) - dL*n(3);
% Test if points are inside or on cone surface
Q    = [(XG - x0),(YG - y0),(ZG - z0)];
F    = Q*n - sqrt(sum(Q.^2,2))*cos(theta);
isin = (F >= -Rfn);
% Subset points
XG  = XG(isin);
YG  = YG(isin);
ZG  = ZG(isin);
IND = IND(isin);
