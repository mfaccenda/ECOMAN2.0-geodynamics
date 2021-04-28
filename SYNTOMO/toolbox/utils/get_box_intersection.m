function [xi,yi,zi] = get_box_intersection(nvec,origin,xmin,xmax,ymin,ymax,zmin,zmax)
% GET_BOX_INTERSECTION: Simple funciton to find intersection of a vector
% with edge of box.
%
% INPUT
%     nvec: Vector components
%   origin: Origin of vector
%     xmin: Minimum x-coordinate of box
%     xmax: Maximum x-coordinate of box
%     ymin: Minimum y-coordinate of box
%     ymax: Maximum y-coordinate of box
%     zmin: Minimum z-coordinate of box
%     zmax: Maximum z-coordinate of box
%
% OUTPUT
%   xi: Intersection x-coordinate
%   yi: Intersection y-coordinate
%   zi: Intersection z-coordinate
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

% Step 1: Towards which quadrant does the vector point?
% Step 2: To which of the three boundaries is the origin closest?
% Step 3: Extend vector towards closest boundary
% Step 4: Check that intersection is not outside box
% Step 5: Try next closest boundary if intersection point is outside box

% Identify which quadrant our vector points to
if (nvec(1) >= 0) && (nvec(2) >= 0) && (nvec(3) >= 0)
    % Upper NE
    bi = [xmax;ymax;zmax] - origin(:);
elseif (nvec(1) < 0) && (nvec(2) >= 0) && (nvec(3) >= 0)
    % Upper NW
    bi = [xmin;ymax;zmax] - origin(:);
elseif (nvec(1) < 0) && (nvec(2) < 0) && (nvec(3) >= 0)
    % Upper SW
    bi = [xmin;ymin;zmax] - origin(:);
elseif (nvec(1) >= 0) && (nvec(2) < 0) && (nvec(3) >= 0)
    % Upper SE
    bi = [xmax;ymin;zmax] - origin(:);
elseif (nvec(1) >= 0) && (nvec(2) >= 0) && (nvec(3) < 0)
    % Lower NE
    bi = [xmax;ymax;zmin] - origin(:);
elseif (nvec(1) < 0) && (nvec(2) >= 0) && (nvec(3) < 0)
    % Lower NW
    bi = [xmin;ymax;zmin] - origin(:);
elseif (nvec(1) < 0) && (nvec(2) < 0) && (nvec(3) < 0)
    % Lower SW
    bi = [xmin;ymin;zmin] - origin(:);
elseif (nvec(1) >= 0) && (nvec(2) < 0) && (nvec(3) < 0)
    % Lower SE
    bi = [xmax;ymin;zmin] - origin(:);
end

% Get intersection
if any(bi == 0)
    xi = [];
    yi = [];
    zi = [];
else
    
    % Initial try
    ibnd = find(abs(bi) == min(abs(bi)),1,'first');
    xi   = origin(1) + bi(ibnd)*nvec(1)/nvec(ibnd);
    yi   = origin(2) + bi(ibnd)*nvec(2)/nvec(ibnd);
    zi   = origin(3) + bi(ibnd)*nvec(3)/nvec(ibnd);
    
    % Check Intersection
    tf_reset = false;
    if (xi < xmin) || (xi > xmax)
        ibnd     = 1;
        tf_reset = true;
    end
    if (yi < ymin) || (yi > ymax)
        ibnd     = 2;
        tf_reset = true;
    end
    if (zi < zmin) || (zi > zmax)
        ibnd     = 3;
        tf_reset = true;
    end
    
    % Recompute intersection point
    if tf_reset
        xi   = origin(1) + bi(ibnd)*nvec(1)/nvec(ibnd);
        yi   = origin(2) + bi(ibnd)*nvec(2)/nvec(ibnd);
        zi   = origin(3) + bi(ibnd)*nvec(3)/nvec(ibnd);
    end
    
    % Final check
    if (xi < xmin) || (xi > xmax) || (yi < ymin) || (yi > ymax) || (zi < zmin) || (zi > zmax)
        error('Bad intersection. Point outside model.');
    end
    
end
