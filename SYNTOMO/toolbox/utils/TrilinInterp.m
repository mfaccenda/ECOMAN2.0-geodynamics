function [g,w,ind] = TrilinInterp(xg,yg,zg,F,xq,yq,zq)
% TRILININTERP: Performs trilinear interpolation weights with option to
% return interpolation weights.
%
% INPUT
%   xg: X-coordinate vector of mesh
%   yg: Y-coordinate vector of mesh
%   zg: Z-ccordinate vector of mesh
%    F: Scalar mesh field
%   xq: X-coordinate(s) of query point
%   yq: Y-coordinate(s) of query point
%   zq: Z-coordinate(s) of query point
%
% OUTPUT
%     g: Value at query point
%     w: Interpolation weight
%   ind: Index of interpolation weight in mesh
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

% Force column vectors
xg = xg(:);
yg = yg(:);
zg = zg(:);

% Get index of cell nearest to query point
[ix,jy,kz] = grid2sub(xg,yg,zg,xq,yq,zq,true);

% Indices of 8 nodes required for trilinear interpolation
i8  = [ix, ix+1, ix+1, ix, ix, ix+1, ix+1, ix];
j8  = [jy, jy, jy+1, jy+1, jy, jy, jy+1, jy+1];
k8  = [kz, kz, kz, kz, kz+1, kz+1, kz+1, kz+1];
ind = sub2ind([length(xg),length(yg),length(zg)],i8,j8,k8);

% Grid coordinates of nodes
x0 = xg(ix);
y0 = yg(jy);
z0 = zg(kz);
x1 = xg(ix+1);
y1 = yg(jy+1);
z1 = zg(kz+1);

% Trilinear interpolation weights
w  = zeros(size(ind));
w(:,1) = (-x1.*y1.*z1 + y1.*z1.*xq + x1.*z1.*yq + x1.*y1.*zq - z1.*xq.*yq...
    - y1.*xq.*zq - x1.*yq.*zq + xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,2) = (x0.*y1.*z1 - y1.*z1.*xq - x0.*z1.*yq - x0.*y1.*zq + z1.*xq.*yq...
    + y1.*xq.*zq + x0.*yq.*zq - xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,3) = (-x0.*y0.*z1 + y0.*z1.*xq + x0.*z1.*yq + x0.*y0.*zq - z1.*xq.*yq...
    - y0.*xq.*zq - x0.*yq.*zq + xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,4) = (x1.*y0.*z1 - y0.*z1.*xq - x1.*z1.*yq - x1.*y0.*zq + z1.*xq.*yq...
    + y0.*xq.*zq + x1.*yq.*zq - xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,5) = (x1.*y1.*z0 - y1.*z0.*xq - x1.*z0.*yq - x1.*y1.*zq + z0.*xq.*yq...
    + y1.*xq.*zq + x1.*yq.*zq - xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,6) = (-x0.*y1.*z0 + y1.*z0.*xq + x0.*z0.*yq + x0.*y1.*zq - z0.*xq.*yq...
    - y1.*xq.*zq - x0.*yq.*zq + xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,7) = (x0.*y0.*z0 - y0.*z0.*xq - x0.*z0.*yq - x0.*y0.*zq + z0.*xq.*yq...
    + y0.*xq.*zq + x0.*yq.*zq - xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));
w(:,8) = (-x1.*y0.*z0 + y0.*z0.*xq + x1.*z0.*yq + x1.*y0.*zq - z0.*xq.*yq...
    - y0.*xq.*zq - x1.*yq.*zq + xq.*yq.*zq)./((x0 - x1).*(y0 - y1).*(z0 - z1));

% Value at query poiont
if ~isempty(F)
    g = sum(w.*F(ind),2);
else
    g = [];
end
