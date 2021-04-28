function K = st_GetHFFKernel(Tc,ray,vray,XG,YG,ZG,IND,tol)
% ST_GETHFFKERNEL: Defines the heuristic finite-frequency travel-time kernel
% (an approximation to the Born sensitivity kernel). Created for SynTomo
% matlab package.
%
% INPUT
%     Tc: Dominant period of wave (s)
%    ray: SynTomo raypath structure (see function 'st_GetRayPath1D')
%   vray: Velocity along ray path
%     XG: Global cartesian x-coordinates at which to compute kernel values.
%     YG: Global cartesian y-coordinates at which to compute kernel values.
%     ZG: Global cartesian z-coordinates at which to compute kernel values.
%    IND: Linear index of XG, YG, and ZG coordinates in model
%
% OUTPUT
%   K: Kernel structure with the following fields...
%         w: Kernel travel-time sensitivity at model points
%       ind: Linear index of weights in model
%       azm: Azimuth of wave propagatioon at kernel points measured CCW
%            from x-axis in RADIANS.
%       elv: Elevation of wave propagatioon at kernel points measured CCW
%            from x,y-plane in RADIANS.
%
% REFERENCES
% + VanderBeek, B.P. & Faccenda, M. (in review). Imaging realistic upper 
%   mantle anisotropy with teleseismic P-wave delays: Insights from 
%   tomographic reconstructions of subduction simulations. Geophysical 
%   Journal International.
%   --> Implementation of approximate finite-frequency kernels for
%       anisotropic models
% + Schmandt, B. & Humphreys, E. (2010). Seismic heterogeneity and 
%   small‐scale convection in the southern California upper mantle. 
%   Geochemistry, Geophysics, Geosystems, 11(5).
%   --> Original description of approximation to Born travel-time
%       sensitivity kernels
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

% Hard-coded Fresnel kernel parameters
% Number of elements defining reference kernel (half-width)
nkern  = 50;
% Kernel smoothing factor (width of gaussian)
ksig   = 0.25;
% Half-width of kernel tapper (number of elements). Taper shape is 
% currently defined as side lobs of unsmoothed kernel.
ktap   = round(nkern*(1-sqrt(1/2)));

% Define the non-dimensional normalized kernel. A scaled version of this
% vector is interpolated to the Fresnel zone slices.
zkern = linspace(-1,1,1 + 2*nkern);
Kw     = sin(pi*(zkern.^2));

% Smooth the kernel
G     = exp(-(1/2)*(zkern/ksig).^2); % Gaussian smoothing
Kw     = conv(G,Kw);
Kw     = Kw((1+nkern):(1+3*nkern));
% Taper kernel. Shape of taper is that describing edges of sin(pi*x^2) function.
W     = [sin(pi*linspace(-1,-sqrt(1/2),ktap).^2),sin(pi*linspace(sqrt(1/2),1,ktap).^2)];
W     = [W(1:ktap),ones(1,1+(2*nkern)-(2*ktap)),W((ktap+1):(2*ktap))];
Kw     = Kw.*W;
% Normalize
Kw     = Kw./max(Kw);

% Reverse along-ray distance such that we measure from the receiver
lray = ray.d(end) - ray.d;

% Define ray unit vector components along path
[n1,n2,n3] = sph2cart(ray.azm,ray.elv,ones(size(ray.azm)));

% Get piecewise linear approximation to ray path
[xpl,ypl,zpl,dpl,ibreak] = break_curve(ray.x,ray.y,ray.z,lray,tol);
xpl = xpl(ibreak);
ypl = ypl(ibreak);
zpl = zpl(ibreak);
% Need to reverse along-ray distance
dpl = dpl(end) - dpl;

% Compute kernel weights
K.w   = zeros(size(XG));
K.azm = zeros(size(XG));
K.elv = zeros(size(XG));
nseg  = length(xpl) - 1;
for ii = 1:nseg
    % Segment unit vectors
    nray_b = [xpl(ii+1)-xpl(ii); ypl(ii+1)-ypl(ii); zpl(ii+1)-zpl(ii)];
    nray_b = nray_b./norm(nray_b);
    if (ii + 1) <= nseg
        nray_a = [xpl(ii+2)-xpl(ii+1); ypl(ii+2)-ypl(ii+1); zpl(ii+2)-zpl(ii+1)];
        nray_a = nray_a./norm(nray_a);
    else
        nray_a = nray_b;
    end
    
    % Select nodes in segment volume
    Fb   = nray_b(1)*(XG - xpl(ii)) + nray_b(2)*(YG - ypl(ii)) + nray_b(3)*(ZG - zpl(ii));
    Fa   = nray_a(1)*(XG - xpl(ii+1)) + nray_a(2)*(YG - ypl(ii+1)) + nray_a(3)*(ZG - zpl(ii+1));
    keep = (Fa <= 0) & (Fb >= 0);
    
    % Compute segment normal distances
    R = sqrt(((XG(keep) - xpl(ii)).^2) + ((YG(keep) - ypl(ii)).^2) + ((ZG(keep) - zpl(ii)).^2) - (Fb(keep).^2));
    
    % Interpolate true along ray distance
    Fb = dpl(ibreak(ii)) - Fb(keep);
    d  = interp1(dpl,lray,Fb,'linear',0);
    v  = interp1(dpl,vray,Fb,'linear',vray(end));
    
    % Interpolate ray orientation
    dx = interp1(dpl,n1,Fb,'linear',1);
    dy = interp1(dpl,n2,Fb,'linear',0);
    dz = interp1(dpl,n3,Fb,'linear',0);
    [phi,theta] = cart2sph(dx,dy,dz);
    
    % Normalized Fresnel radius
    Rf = sqrt(Tc*v.*d.*(ray.L - d)./ray.L)/sqrt(2);
    R  = R./Rf;
    
    % Kernel weights
    K.w(keep)   = interp1(zkern,Kw,R,'linear',0);
    K.w(keep)   = K.w(keep)./(pi*(Rf.^2));
    K.azm(keep) = phi;
    K.elv(keep) = theta;
end
% Re-normalize kernel weights
K.w   = lray(1)*K.w./sum(K.w(:));
K.ind = IND;
% Only need to return non-zero weights
irmv        = (K.w == 0);
K.w(irmv)   = [];
K.ind(irmv) = [];
K.azm(irmv) = [];
K.elv(irmv) = [];
