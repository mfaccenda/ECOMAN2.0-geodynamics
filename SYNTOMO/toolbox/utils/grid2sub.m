function [ix,jy,kz] = grid2sub(xg,yg,zg,xq,yq,zq,tf_cell)
% GRID2SUB: Convert cartesian coordinates to nearest grid point subscript 
% in regular mesh.
%
% INPUT
%        xg: X-coordinate grid vector
%        yg: Y-coordinate grid vector
%        zg: Z-coordinate grid vector
%        xq: X-coordinate of query point(s)
%        yq: Y-coordinate of query point(s)
%        zq: Z-coordinate of query point(s)
%   tf_cell: If true, returns subscript of nearest cell rather than node
%
% OUPUT
%   ix: First dimension index
%   jy: Second dimension index
%   kz: Third dimension index
%
% B. VanderBeek
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

% Convert mesh coordinate vectors to cell coordinate vectors?
if tf_cell
    xg = (diff(xg)./2) + xg(1:end-1);
    yg = (diff(yg)./2) + yg(1:end-1);
    zg = (diff(zg)./2) + zg(1:end-1);
end

% Find nearest grid point in mesh
ix = zeros(size(xq));
jy = zeros(size(yq));
kz = zeros(size(zq));
% This loop could be vectorized to improve speed.
for n = 1:length(xq)
    dx = abs(xg - xq(n));
    ix(n) = find(dx == min(dx),1,'first');
    
    dy = abs(yg - yq(n));
    jy(n) = find(dy == min(dy),1,'first');
    
    dz = abs(zg - zq(n));
    kz(n) = find(dz == min(dz),1,'first');
end
