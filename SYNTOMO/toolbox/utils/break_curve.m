function [xpl,ypl,zpl,dpl,ibreak] = break_curve(xcrv,ycrv,zcrv,lcrv,tol)
% BREAK_CURVE: Breaks a curve into piece-wise linear segments. Created for 
% SynTomo matlab package.
%
% INPUT
%   xcrv: Vector defining curve x-coordinates
%   ycrv: Vector defining curve y-coordinates
%   zcrv: Vector defining curve z-coordinates
%   lcrv: Distance along curve
%    tol: Maximum allowed distance between linear segment and true curve
%
% OUTPUT
%      xpl: X-coordinates of piece-wise linear curve
%      ypl: Y-coordinates of piece-wise linear curve
%      zpl: Z-coordinates of piece-wise linear curve
%      dpl: Distance along piece-wise linear curve
%   ibreak: Starting indices of linear segments
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

% Initial approximation to ray path--straight line
npts = length(xcrv);
% Slopes
dxdl = (xcrv(end) - xcrv(1))/(lcrv(end) - lcrv(1));
dydl = (ycrv(end) - ycrv(1))/(lcrv(end) - lcrv(1));
dzdl = (zcrv(end) - zcrv(1))/(lcrv(end) - lcrv(1));
% Intercepts
x0 = xcrv(1) - dxdl*lcrv(1);
y0 = ycrv(1) - dydl*lcrv(1);
z0 = zcrv(1) - dzdl*lcrv(1);
% Linear ray path
xpl = dxdl*lcrv + x0;
ypl = dydl*lcrv + y0;
zpl = dzdl*lcrv + z0;
% Error in approximation
ep  = sqrt(((xpl - xcrv).^2) + ((ypl - ycrv).^2) + ((zpl - zcrv).^2));
ib1 = 1;
ib2 = find(ep == max(ep),1,'first');
ib3 = npts;
ep  = max(ep);

% Track segments
iseg = ones(npts,1);
nseg = 1;
% Break the ray path until tolerence is met by all segments
while (ep > tol) && (nseg < round(npts/2))
    % First new segment
    dxdl = (xcrv(ib2) - xcrv(ib1))/(lcrv(ib2) - lcrv(ib1));
    dydl = (ycrv(ib2) - ycrv(ib1))/(lcrv(ib2) - lcrv(ib1));
    dzdl = (zcrv(ib2) - zcrv(ib1))/(lcrv(ib2) - lcrv(ib1));
    
    x0 = xcrv(ib1) - dxdl*lcrv(ib1);
    y0 = ycrv(ib1) - dydl*lcrv(ib1);
    z0 = zcrv(ib1) - dzdl*lcrv(ib1);
    
    xpl(ib1:ib2)  = dxdl*lcrv(ib1:ib2) + x0;
    ypl(ib1:ib2)  = dydl*lcrv(ib1:ib2) + y0;
    zpl(ib1:ib2)  = dzdl*lcrv(ib1:ib2) + z0;
    iseg(ib1:ib2) = nseg;
    
    % Second new segment
    dxdl = (xcrv(ib3) - xcrv(ib2))/(lcrv(ib3) - lcrv(ib2));
    dydl = (ycrv(ib3) - ycrv(ib2))/(lcrv(ib3) - lcrv(ib2));
    dzdl = (zcrv(ib3) - zcrv(ib2))/(lcrv(ib3) - lcrv(ib2));
    
    x0 = xcrv(ib2) - dxdl*lcrv(ib2);
    y0 = ycrv(ib2) - dydl*lcrv(ib2);
    z0 = zcrv(ib2) - dzdl*lcrv(ib2);
    
    xpl(ib2+1:ib3)  = dxdl*lcrv(ib2+1:ib3) + x0;
    ypl(ib2+1:ib3)  = dydl*lcrv(ib2+1:ib3) + y0;
    zpl(ib2+1:ib3)  = dzdl*lcrv(ib2+1:ib3) + z0;
    iseg(ib2+1:ib3) = nseg + 1;
    
    % New error and break point
    ep  = sqrt(((xpl - xcrv).^2) + ((ypl - ycrv).^2) + ((zpl - zcrv).^2));
    ib2 = find(ep == max(ep),1,'first');
    ep  = max(ep);
    ib1 = find(iseg == iseg(ib2),1,'first');
    ib3 = find(iseg == iseg(ib2),1,'last');
    
    % Update segment counter
    nseg = nseg + 1;
end
% Check that we did not reach segment limit
if (nseg >= round(npts/2))
    warning('Maximum number of segments reached in ''break_curve''.');
end
% Distance along piecewise linear curve
dpl = cumsum(sqrt((diff(xpl).^2) + (diff(ypl).^2) + (diff(zpl).^2)));
dpl = cat(1,0,dpl);

% Segment indices
ibreak = find(abs(diff(iseg)) > 0);
ibreak = cat(1,1,ibreak+1,npts);

