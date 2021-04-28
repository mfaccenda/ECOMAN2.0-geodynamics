function C = get_cross_section(xg,yg,zg,F,slice)
% GET_CROSS_SECTION: Interpolates 3D model onto a user-specified
% horizontally- or vertically-oriented cross-section plane. Created for 
% SynTomo matlab package.
%
% INPUT
%      xg: x1-coordinate vector of model array
%      yg: x2-coordinate vector of model array
%      zg: x3-coordinate vector of model array
%       F: Model array
%   slice: For a horizontal cross-section plane, 'slice' specifies depth of
%          cross-section. For a vertically oriented cross-section plane,
%          'slice' is a 1x3 vector defining the x,y-origin and azimuth
%          (degrees CCW of x-axis) of the cross-section plane.
%
% OUTPUT
%   C: Cross-section structure with the following fields...
%       F: Model values on cross-section plane
%      X1: Cross-section x-coordinate array
%      X2: Cross-section y-coordinate array
%       X: Model x-coordinates of cross-section
%       Y: Model y-coordinates of cross-section
%       Z: Model z-coordinates of cross-section
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

% Depth slice or cross-section parameters
if length(slice) == 1
    tf_depth = true;
elseif length(slice) == 3
    tf_depth = false;
    a  = slice(1);
    b  = slice(2);
    az = slice(3);
else
    error(['The input variable ''slice'' must either be a scalar for a ',...
        'depth slice or a vector of length 3 for a cross section.']);
end

if tf_depth
    % Interpolate field to desired depth
    F = interp3(yg(:)',xg(:),zg(:)',F,yg(:)',xg(:),slice);
    
    % Define model coordinates in cross-section
    [Y,X] = meshgrid(yg,xg);
    Z     = slice*ones(size(X));
    
    % Coordinate arrays for cross-section plotting
    X1 = X;
    X2 = Y;
else
    % Cross-section
    % Get intersection of cross-section with model boundaries
    [dx1,dy1,dx2,dy2] = get_cross_intersect(min(xg),max(xg),min(yg),max(yg),a,b,az);
    
    % Define model coordinates in cross-section
    dr = (mean(diff(xg)) + mean(diff(yg)))/2; % Horizontal spacing
    L  = distance(dy1,dx1,dy2,dx2);
    n  = 1 + round(L/dr); % Number of points in x1-direction
    X  = linspace(a+dx2,a+dx1,n); % True x-coordinates of cross-section
    Y  = linspace(b+dy2,b+dy1,n); % True y-coordinates of cross-section
    
    % Interpolate model to cross-section plane
    X = repmat(X(:),1,length(zg));
    Y = repmat(Y(:),1,length(zg));
    Z = repmat(zg(:)',n,1);
    F = interp3(yg(:)',xg(:),zg(:)',F,Y,X,Z);
    
    % Coordinate arrays for cross-section
    X1 = distance(b,a,Y,X);
    X2 = Z;
    
    % Define negative coordinates
    ineg = find(X1(:,1) == min(X1(:,1)),1,'first');
    X1(1:ineg,:) = -X1(1:ineg,:);
    
    % Map to global cartesian coordinates
    [X1,~,X2] = geo2cart(X1-mean(X1(:)),zeros(size(X1)),X2 + 6371,0,0);
end

% Ouput structure
C.F   = F;
C.X1  = X1;
C.X2  = X2;
C.X   = X;
C.Y   = Y;
C.Z   = Z;
end

%% Function for identifying intersection points

function [dx1,dy1,dx2,dy2] = get_cross_intersect(xmin,xmax,ymin,ymax,a,b,az)
% Check that slice origin is in the model
if (a > xmax) || (a < xmin)
    error('Angular slice x-origin is outside model domain');
end
if (b > ymax) || (b < ymin)
    error('Angular slice y-origin is outside model domain');
end

% Identify intersection with model boundaries
if sign(tand(az)) >= 0
    dy1j = (ymax - b);
    dy2j = (ymin - b);
else
    dy1j = (ymin - b);
    dy2j = (ymax - b);
end
% Intersection 1 possibilities
dx1i = (xmax - a);
dy1i = dx1i*tand(az);
dx1j = dy1j/tand(az);
% Intersection 2 possibilities
dx2i = (xmin - a);
dy2i = dx2i*tand(az);
dx2j = dy2j/tand(az);
% Locate intersections--point 1
if ((b + dy1i) <= ymax) && ((b + dy1i) >= ymin)
    dx1 = dx1i;
    dy1 = dy1i;
elseif ((a + dx1j) <= xmax) && ((a + dx1j) >= xmin)
    dx1 = dx1j;
    dy1 = dy1j;
end
% Locate intersections--point 2
if ((b + dy2i) <= ymax) && ((b + dy2i) >= ymin)
    dx2 = dx2i;
    dy2 = dy2i;
elseif ((a + dx2j) <= xmax) && ((a + dx2j) >= xmin)
    dx2 = dx2j;
    dy2 = dy2j;
end
end