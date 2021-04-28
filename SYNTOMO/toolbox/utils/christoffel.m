function [Vp,Vss,Vsf,pP,pSs,pSf,X,Y,Z] = christoffel(Cij,rho,numpts,varargin)
% CHRISTOFFEL: Calculate directional dependence of seismic wave speeds and
% polarizations using Christoffel equations.
%
% INPUT
%         Cij: Voigt representation of elastic tensor in GPa (6x6 matrix)
%         rho: Density in kg/m^3 (scalar)
%      numpts: Three options. (1) A scalar specifying the number of points 
%              defining a unit sphere at which wave speeds and
%              polarizations will be computed. (2) A n-by-2 array 
%              specifying the azimuth and elevation of propagation. (3) A
%              n-by-3 array specifying the x, y, and z components of the
%              unit vector pointing in the direction of propagation.
%              -- Azimuth is measured CCW from x1-axis; Elevation is
%                 measured CCW from the x,y-plane. Angles are in DEGREES.
% <Optional>
%   tf_PlotOn: If true, function makes velocity surface plot if numpts is a
%              scalar (option 1); default is true when numpts is a scalar.
%
% OUTPUT
%    Vp: Compressional velocity in direction of propagation (km/s)
%   Vss: Shear velocity in slow polarization direction (km/s)
%   Vsf: Shear velocity in fast polarization direction (km/s)
%    pP: P-wave polarization vector
%   pSs: Slow S-wave polarization vector
%   pSf: Fast S-wave polarization vector
%     X: X-component of propagation direction unit vector
%     Y: Y-component of propagation direction unit vector
%     Z: Z-component of propagation direction unit vector
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


% Variable input arguments
tf_PlotOn = true;
if length(varargin) == 1
    tf_PlotOn = varargin{1};
elseif length(varargin) > 1
    error('Too many optional input arguments.');
end
% Only make plot if input points define a sphere
if ~isscalar(numpts)
    tf_PlotOn = false;
end

% Elastic tensor in Pa
Cij = Cij*(10^9);

% Define propagation directions at which to calculate velocities
if isscalar(numpts)
    % Unit sphere
    [X,Y,Z] = sphere(numpts);
    nout    = (numpts+1)^2;
    % Vectorize
    X = X(:);
    Y = Y(:);
    Z = Z(:);
elseif size(numpts,2) == 2
    % User-specified angular directions
    numpts  = numpts.*pi./180;
    [X,Y,Z] = sph2cart(numpts(:,1),numpts(:,2),ones(size(numpts(:,1))));
    nout    = length(X);
elseif size(numpts,2) == 3
    % User-specified vectoral directions
    X    = numpts(:,1);
    Y    = numpts(:,2);
    Z    = numpts(:,3);
    nout = length(X);
end

% Pre-allocate result arrays
Vp  = zeros(nout,1);
Vss = zeros(nout,1);
Vsf = zeros(nout,1);
pP  = zeros(nout,3);
pSs = zeros(nout,3);
pSf = zeros(nout,3);

% Loop over propagation directions
for m = 1:nout
    % Form Cristoffel matrix without performing 4-fold loop. This 
    % definition was taken from the function 'MS_phasevels' provided as 
    % part of the MSAT (Matlab Seismic Anisotropy Toolbox) package
    % (Walker & Wookey, Computers & Geoscience 2012). Original reference
    % for this equation is:
    % + Winterstein D. F. (1990). Velocity anisotropy terminology for 
    %   geophysicists. Geophysics, vol 55, pp1070-1088.
    N = [X(m) 0     0    0    Z(m) Y(m);...
         0    Y(m)  0    Z(m) 0    X(m);...
         0    0     Z(m) Y(m) X(m) 0];
    G = N*(Cij./rho)*N';
    
    % Solve Eigenvalue problem for wave speeds and polarizations
    [V,L] = eig(G,'vector');
    L     = (1/1000)*sqrt(L(:));
    
    % Define phase velocities
    % P-velocity must be the largest
    ip = (L == max(L));
    % S-slow velocity must be the minimum
    iss = (L == min(L));
    % This handles situation where only one S-wave exists
    if (sum(iss) == 2)
        % If shear wave speeds are the same, then just default to the
        % first entry
        iss(find(iss,1,'last')) = false;
    elseif (sum(iss) == 3) || (sum(iss) == 0)
        error('Problem identifying shear phase.');
    end
    % S-fast must be...
    isf = ~ip & ~iss;
    
    % Check for unique phase identifications
    if ~((sum(ip) == 1) && (sum(iss) == 1) && (sum(isf) == 1))
        error('Problem making phase identification.');
    end
    
    % Use consistent sign convention. There is a +/- 180 ambiguity in the
    % polarization direction. Here we are defining the eigen vectors such
    % that the P-wave polarization orientation is nearest to the
    % propagation direction.
    pn = dot(V(:,ip),N(:,1));
    pr = dot(-V(:,ip),N(:,1));
    if pr > pn
        V = -V;
    end
    
    % Store results
    Vp(m)    = L(ip);
    Vss(m)   = L(iss);
    Vsf(m)   = L(isf);
    pP(m,:)  = V(:,ip)';
    pSs(m,:) = V(:,iss)';
    pSf(m,:) = V(:,isf)';
end
% Reshape if constructing velocity sphere
if isscalar(numpts)
    X   = reshape(X,numpts+1,numpts+1);
    Y   = reshape(Y,numpts+1,numpts+1);
    Z   = reshape(Z,numpts+1,numpts+1);
    Vp  = reshape(Vp,numpts+1,numpts+1);
    Vss = reshape(Vss,numpts+1,numpts+1);
    Vsf = reshape(Vsf,numpts+1,numpts+1);
end

% Plots velocity surfaces
if tf_PlotOn
    
    % P-wave surface
    figure;
    surf(Vp.*X,Vp.*Y,Vp.*Z,Vp,'edgecolor','none'); axis image; colorbar
    caxis([min(Vp(:)),max(Vp(:))])
    xlabel('V_1 (km/s)'); ylabel('V_2 (km/s)'); zlabel('V_3 (km/s)');
    title('P-wave Velocity Surface');
    
    % S-slow surface
    figure;
    surf(Vss.*X,Vss.*Y,Vss.*Z,Vss,'edgecolor','none'); axis image; colorbar
    caxis([min(Vss(:)),max(Vss(:))])
    xlabel('V_1 (km/s)'); ylabel('V_2 (km/s)'); zlabel('V_3 (km/s)');
    title('S-slow Velocity Surface');
    
    % S-fast surface
    figure;
    surf(Vsf.*X,Vsf.*Y,Vsf.*Z,Vsf,'edgecolor','none'); axis image; colorbar
    caxis([min(Vsf(:)),max(Vsf(:))])
    xlabel('V_1 (km/s)'); ylabel('V_2 (km/s)'); zlabel('V_3 (km/s)');
    title('S-fast Velocity Surface');
    
end
 