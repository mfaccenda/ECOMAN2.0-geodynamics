function [L,is,js] = LaplaceCoeff(x,ix,Nijk)
% LAPLACECOEFF: Defines finite-difference Laplacian coefficients for 1D,
% 2D, or 3D arrays.
%
% INPUT
%      x: A grid vector.
%     ix: Specifies which dimension the grid vector 'x' is defining
%   Nijk: A 3x1 array defining the number of elements in each dimension
%
% OUTPUT
%    L: Laplacian coefficients
%   is: Row index in array
%   js: Column index in array
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

% Number of elements in input dimension
nx = length(x);

% Number of block diagonal elements required
m = Nijk(1)*Nijk(2)*Nijk(3)/nx;

if nx > 2
    % Forward Laplacian coefficients
    f   = (x(1) - 2*x(2) + x(3))/(x(2)-x(1));
    wif = [1,-(2+f)/(1+f),1/(1+f)]; % ii, ii+1, ii+2
    
    % Central Laplacian coefficients
    dxf = x(3:end) - x(2:end-1);
    dxb = x(2:end-1) - x(1:end-2);
    f   = (dxf - dxb)./(dxf + dxb);
    wic = [-(1+f(:))./2,ones(nx-2,1),-(1-f(:))./2]; % ii-1, ii, ii+1
    
    % Backward Laplacian coefficients
    f   = (x(end) - 2*x(end-1) + x(end-2))/(x(end)-x(end-1));
    wib = [1/(1+f),-(2+f)/(1+f),1]; % ii-2, ii-1, ii
    
    % Build coefficient and indexing arrays for the first nx-points
    L  = cat(1,wif,wic,wib);
    % The ith-subscripts for the coefficients
    ri = repmat((1:nx),3,1)';
    ci = ri + cat(1,[0,1,2],repmat([-1,0,1],nx-2,1),[-2,-1,0]);
    
    % Replicate for each point in the remaining dimensions
    L  = repmat(L,m,1);
    ri = repmat(ri,m,1);
    ci = repmat(ci,m,1);
    
    % Define subscripts in the other two dimensions
    nsub      = Nijk;
    nsub(ix)  = [];
    [rj,~,rk] = meshgrid(1:nsub(1),(1:nx),1:nsub(2));
    rj        = rj(:);
    rk        = rk(:);
    
    % Define sparse matrix subscripts and weights
    if ix == 1
        is = sub2ind(Nijk,ri,repmat(rj,1,3),repmat(rk,1,3));
        js = sub2ind(Nijk,ci,repmat(rj,1,3),repmat(rk,1,3));
    elseif ix == 2
        is = sub2ind(Nijk,repmat(rj,1,3),ri,repmat(rk,1,3));
        js = sub2ind(Nijk,repmat(rj,1,3),ci,repmat(rk,1,3));
    elseif ix == 3
        is = sub2ind(Nijk,repmat(rj,1,3),repmat(rk,1,3),ri);
        js = sub2ind(Nijk,repmat(rj,1,3),repmat(rk,1,3),ci);
    end
    
    % Return column vector
    is = is(:);
    js = js(:);
    L  = L(:);
else
    L  = [];
    is = [];
    js = [];
end
