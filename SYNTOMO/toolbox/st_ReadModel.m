function model = st_ReadModel(theFile,nq1,nq2,nq3,tf_cart,varargin)
% ST_READMODEL: Reads in elastic model created via VIZTOMO, a part of the 
% ECOMAN software package, and interpolates to user-specified grid. Output
% is a SynTomo model structure.
%
% INPUT
%   theFile: Name of VIZTOMO HDF5 model file to be read.
%       nq1: Desired number of nodes in 1st (longitudinal/x) dimension of
%            SynTomo model.
%       nq2: Desired number of nodes in 2nd (latitudinal/y) dimension of
%            SynTomo model.
%       nq3: Desired number of nodes in 3rd (radial/z) dimension of
%            SynTomo model.
%   tf_cart: If true, the VIZTOMO model file is in cartesian coordinates.
%            If false, the VIZTOMO model is in spherical coordinates. See
%            below for assumptions regarding each coordinate system.
% <optional>
%     order: A 1x3 vector specifying the array dimensions corresponding to
%            logitudinal/x-, latitudinal/y-, and radial/z-directions. 
%            Default order is assumed to be [1,2,3].
%
% OUTPUT
%   model: SynTomo model structure with the following fields...
%            long: Longitude coordinate vector (deg.)
%            latg: Latitude coordinate vector (deg.)
%             rzg: Depth coordinate vector; negative into Earth with 
%                  surface at 0 (km).
%            wlon: Width of model in longitudinal direction (deg.)
%            wlat: Width of model in latitudinal direction (deg.)
%           minrz: Minimum depth of array (km); Maximum depth corresponds
%                  to Earth's surface at rzg = 0 km.
%              n1: Number of nodes in first dimension
%              n2: Number of nodes in second dimension
%              n3: Number of nodes in third dimension
%             RHO: Density model (kg/m^3)
%             Cij: One field for each of the 21 elastic coefficients (GPa)
%              Vp: Average compressional wave speed (km/s)
%              Vs: Average shear wave speed (km/s)
%              Re: Spherical model radius (km)
%
% VIZTOMO COORDINATE SYSTEM
% + Spherical models: The 1st, 2nd, and 3rd model dimensions correspond to
%   longitude (radians), latitide (radians), and radius (m), respectively,
%   unless otherwise specified by the optional input 'order'.
% + Cartesian models: The 1st, 2nd, and 3rd model dimensions correspond to
%   x (m), y (m), and z (m), respectively, unless otherwise specified by 
%   the optional input 'order'. Additionally, z is assumed to define depth
%   positive into the Earth.
% + Reference Earth radius is assumed to be 6371 km.
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

% Define optional inputs
if isempty(varargin)
    order = [1,2,3];
elseif length(varargin) == 1
    order = varargin{1};
else
    error('Too many optional input arguments');
end

% Read grid vectors
x1 = h5read(theFile,['/X',num2str(order(1))]);
x2 = h5read(theFile,['/X',num2str(order(2))]);
x3 = h5read(theFile,['/X',num2str(order(3))]);
% Original dimensions
n  = [length(x1),length(x2),length(x3)];
n  = n(order);

% Desired grid vectors
xq1 = linspace(x1(1),x1(end),nq1)';
xq2 = linspace(x2(1),x2(end),nq2)';
xq3 = linspace(x3(1),x3(end),nq3)';

% Store grid information in model structure. Origin of grid is made to
% correspond to center of grid.
model.long   = xq1 - ((max(xq1) + min(xq1))/2);
model.latg   = xq2 - ((max(xq2) + min(xq2))/2);
model.rzg    = xq3;
model.wlon   = max(model.long) - min(model.long);
model.wlat   = max(model.latg) - min(model.latg);
model.minrz  = min(model.rzg);
model.n1     = nq1;
model.n2     = nq2;
model.n3     = nq3;

% Read in density model
model.RHO = permute(reshape(h5read(theFile,'/Rho'),n(1),n(2),n(3)),order);

% Read in elasticity model
Cij = h5read(theFile,'/Sav');

% Re-arrange coefficient order
row = [1,6,5,10,14,18,2,4,16,11,15,3,13,17,12,7,21,20,8,19,9];
Cij = Cij(row,:);

% Reshape and store elastic components
model.C11 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C12 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C13 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C14 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C15 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C16 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C22 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C23 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C24 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C25 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C26 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C33 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C34 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C35 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C36 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C44 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C45 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C46 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C55 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C56 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); Cij(1,:) = [];
model.C66 = permute(reshape(Cij(1,:),n(1),n(2),n(3)),order); clear('Cij');

% Need to rotate tensors for coordinate system changes
if all(order == [1,2,3])
    % Nothing to do in this case but considered for completness-sake
elseif all(order == [1,3,2])
    % A 90-degree rotation about x-axis
    % C22/C33
    Cij = model.C22;
    model.C22 = model.C33;
    model.C33 = Cij;
    % C55/C66
    Cij = model.C55;
    model.C55 = model.C66;
    model.C66 = Cij;
    % C12/C13
    Cij = model.C12;
    model.C12 = model.C13;
    model.C13 = Cij;
    % C14
    model.C14 = -model.C14;
    % C15/C16
    Cij = model.C15;
    model.C15 = -model.C16;
    model.C16 = Cij;
    % C24/C34
    Cij = model.C24;
    model.C24 = -model.C34;
    model.C34 = -Cij;
    % C25/C36
    Cij = model.C25;
    model.C25 = -model.C36;
    model.C36 = Cij;
    % C26/C35
    Cij = model.C26;
    model.C26 = model.C35;
    model.C35 = -Cij;
    % C45/C46
    Cij = model.C45;
    model.C45 = model.C46;
    model.C46 = -Cij;
    % C56
    model.C56 = -model.C56;
elseif all(order == [2,1,3])
    error('Rotation not implemented.');
elseif all(order == [2,3,1])
    error('Rotation not implemented.');
elseif all(order == [3,1,2])
    error('Rotation not implemented.');
elseif all(order == [3,2,1])
    error('Rotation not implemented.');
else
    error('Undefined case.');
end

% Need to convert model to spherical coordinates?
if tf_cart
    x3 = -x3; % SynTomo models define z negative into Earth
    [x1,x2,model,xq1,xq2,xq3] = st_ModelCart2Geo(x1,x2,x3,n,order,model,nq1,nq2,nq3);
end

% Check if we need to interpolate
if (model.n1 == n(order(1))) && (model.n2 == n(order(2))) && (model.n3 == n(order(3)))
    tf_interp = (max(abs(xq1 - x1)) > 10*eps) || (max(abs(xq2 - x2)) > 10*eps) || (max(abs(xq3 - x3)) > 10*eps);
else
    tf_interp = true;
end

% Interpolate model to desired grid spacing
if tf_interp
    fldnm = fieldnames(model);
    for ifld = 1:length(fldnm)
        if (size(model.(fldnm{ifld}),1) == n(order(1))) &&...
                (size(model.(fldnm{ifld}),2) == n(order(2))) &&...
                (size(model.(fldnm{ifld}),3) == n(order(3)))
            model.(fldnm{ifld}) = interp3(x2(:)',x1(:),x3(:)',...
                model.(fldnm{ifld}),xq2(:)',xq1(:),xq3(:)');
            if any(isnan(model.(fldnm{ifld})(:)))
                error('NaN values after model interpolation.');
            end
        end
    end
end

% % Interpolated elastic coefficients
% if tf_interp
%     model.RHO = interp3(x2(:)',x1(:),x3(:)',model.RHO,xq2(:)',xq1(:),xq3(:)');
%     model.C11 = interp3(x2(:)',x1(:),x3(:)',model.C11,xq2(:)',xq1(:),xq3(:)');
%     model.C12 = interp3(x2(:)',x1(:),x3(:)',model.C12,xq2(:)',xq1(:),xq3(:)');
%     model.C13 = interp3(x2(:)',x1(:),x3(:)',model.C13,xq2(:)',xq1(:),xq3(:)');
%     model.C14 = interp3(x2(:)',x1(:),x3(:)',model.C14,xq2(:)',xq1(:),xq3(:)');
%     model.C15 = interp3(x2(:)',x1(:),x3(:)',model.C15,xq2(:)',xq1(:),xq3(:)');
%     model.C16 = interp3(x2(:)',x1(:),x3(:)',model.C16,xq2(:)',xq1(:),xq3(:)');
%     model.C22 = interp3(x2(:)',x1(:),x3(:)',model.C22,xq2(:)',xq1(:),xq3(:)');
%     model.C23 = interp3(x2(:)',x1(:),x3(:)',model.C23,xq2(:)',xq1(:),xq3(:)');
%     model.C24 = interp3(x2(:)',x1(:),x3(:)',model.C24,xq2(:)',xq1(:),xq3(:)');
%     model.C25 = interp3(x2(:)',x1(:),x3(:)',model.C25,xq2(:)',xq1(:),xq3(:)');
%     model.C26 = interp3(x2(:)',x1(:),x3(:)',model.C26,xq2(:)',xq1(:),xq3(:)');
%     model.C33 = interp3(x2(:)',x1(:),x3(:)',model.C33,xq2(:)',xq1(:),xq3(:)');
%     model.C34 = interp3(x2(:)',x1(:),x3(:)',model.C34,xq2(:)',xq1(:),xq3(:)');
%     model.C35 = interp3(x2(:)',x1(:),x3(:)',model.C35,xq2(:)',xq1(:),xq3(:)');
%     model.C36 = interp3(x2(:)',x1(:),x3(:)',model.C36,xq2(:)',xq1(:),xq3(:)');
%     model.C44 = interp3(x2(:)',x1(:),x3(:)',model.C44,xq2(:)',xq1(:),xq3(:)');
%     model.C45 = interp3(x2(:)',x1(:),x3(:)',model.C45,xq2(:)',xq1(:),xq3(:)');
%     model.C46 = interp3(x2(:)',x1(:),x3(:)',model.C46,xq2(:)',xq1(:),xq3(:)');
%     model.C55 = interp3(x2(:)',x1(:),x3(:)',model.C55,xq2(:)',xq1(:),xq3(:)');
%     model.C56 = interp3(x2(:)',x1(:),x3(:)',model.C56,xq2(:)',xq1(:),xq3(:)');
%     model.C66 = interp3(x2(:)',x1(:),x3(:)',model.C66,xq2(:)',xq1(:),xq3(:)');
% end

% Flip model z-coordinate
if x3(end) > x3(1)
    fldnm = fieldnames(model);
    for ifld = 1:length(fldnm)
        if (size(model.(fldnm{ifld}),1) == model.n1) &&...
                (size(model.(fldnm{ifld}),2) == model.n2) &&...
                (size(model.(fldnm{ifld}),3) == model.n3)
            model.(fldnm{ifld}) = model.(fldnm{ifld})(:,:,end:-1:1);
        end
    end
    model.rzg = model.rzg(end:-1:1);
end

% Isotropic Velocities
model.Vp = (model.C11 + model.C22 + model.C33 + 2*model.C12 + 2*model.C13 + 2*model.C23)./9;
model.Vs = (2*model.C11 + 2*model.C22 + 2*model.C33 + 6*model.C44 + 6*model.C55...
    + 6*model.C66 - 2*model.C12 - 2*model.C13 - 2*model.C23)./30;
model.Vp = sqrt((10^9)*(model.Vp + (4/3)*model.Vs)./model.RHO)./1000;
model.Vs = sqrt((10^9)*model.Vs./model.RHO)./1000;

% Convert coordinates to degrees/kilometers
model.Re     = max(model.rzg(:))/1000;
model.long   = 180*model.long./pi;
model.latg   = 180*model.latg./pi;
model.rzg    = (model.rzg./1000) - model.Re;
model.wlon   = 180*model.wlon/pi;
model.wlat   = 180*model.wlat/pi;
model.minrz  = (model.minrz/1000) - model.Re;
end

%% Convert Cartesian Model to Spherical
function [x1,x2,model,xq1,xq2,xq3] = st_ModelCart2Geo(x1,x2,x3,n,order,model,nq1,nq2,nq3)

% Earth radius
Re = 6371;

% Map cartesian model coordinates to geographic coordinates assuming model
% is centered at (0N,0E).
xkm       = (x1 - ((max(x1) + min(x1))/2))./1000;
ykm       = (x2 - ((max(x2) + min(x2))/2))./1000;
elv       = x3./1000;
[Y,X]     = meshgrid(ykm,xkm);
[LON,LAT] = cart2geo(X,Y,zeros(size(X)),0,0);

% Need to interpolate model to regular array in geographic coordinates.
% Define the cartesian coordinates of this regular geographic array for
% interpolation.
wlon = min(abs(LON(1,:)));
wlat = min(abs(LAT(:,1)));
long = linspace(-wlon,wlon,n(order(1)));
latg = linspace(-wlat,wlat,n(order(2)));
[Yq,Xq,Zq]   = meshgrid(latg,long,elv);
[Xq,Yq,~,Zq] = geo2cart(Xq,Yq,Zq + Re,0,0); % Here, Zq is elevation
% Fix rounding errors to prevent NaN's in interpolation
Zq = min(Zq,0);
Zq = max(Zq,min(elv));

% Interpolate model to regular geographic grid
fldnm = fieldnames(model);
for ifld = 1:length(fldnm)
    if (size(model.(fldnm{ifld}),1) == n(order(1))) &&...
            (size(model.(fldnm{ifld}),2) == n(order(2))) &&...
            (size(model.(fldnm{ifld}),3) == n(order(3)))
        model.(fldnm{ifld}) = interp3(ykm(:)',xkm(:),elv(:)',model.(fldnm{ifld}),Yq,Xq,Zq);
        if any(isnan(model.(fldnm{ifld})(:)))
            error('NaN values after model interpolation to regular geographic grid.');
        end
    end
end

% Rotate tensors for spherical geometry
for m = 1:prod(n)
    % The mth Voigt tensor
    C = [model.C11(m),model.C12(m),model.C13(m),model.C14(m),model.C15(m),model.C16(m);...
        model.C12(m),model.C22(m),model.C23(m),model.C24(m),model.C25(m),model.C26(m);...
        model.C13(m),model.C23(m),model.C33(m),model.C34(m),model.C35(m),model.C36(m);...
        model.C14(m),model.C24(m),model.C34(m),model.C44(m),model.C45(m),model.C46(m);...
        model.C15(m),model.C25(m),model.C35(m),model.C45(m),model.C55(m),model.C56(m);...
        model.C16(m),model.C26(m),model.C36(m),model.C46(m),model.C56(m),model.C66(m)];
    [ii,jj,~] = ind2sub(n(order),m);
    C = rotate_tensor_voigt(C,latg(jj),-long(ii),0,[1,2,3]);
    
    % Update elastic model
    % Row 1
    model.C11(m) = C(1,1);
    model.C12(m) = C(1,2);
    model.C13(m) = C(1,3);
    model.C14(m) = C(1,4);
    model.C15(m) = C(1,5);
    model.C16(m) = C(1,6);
    % Row 2
    model.C22(m) = C(2,2);
    model.C23(m) = C(2,3);
    model.C24(m) = C(2,4);
    model.C25(m) = C(2,5);
    model.C26(m) = C(2,6);
    % Row 3
    model.C33(m) = C(3,3);
    model.C34(m) = C(3,4);
    model.C35(m) = C(3,5);
    model.C36(m) = C(3,6);
    % Row 4
    model.C44(m) = C(4,4);
    model.C45(m) = C(4,5);
    model.C46(m) = C(4,6);
    % Row 5
    model.C55(m) = C(5,5);
    model.C56(m) = C(5,6);
    % Row 6
    model.C66(m) = C(6,6);
end

% Update parameters to reflect coordinate system change
% Initial grid vectors
x1 = long*pi./180;
x2 = latg*pi./180;
% Desired grid vectors
xq1 = linspace(x1(1),x1(end),nq1)';
xq2 = linspace(x2(1),x2(end),nq2)';
xq3 = linspace(x3(1),x3(end),nq3)';
% Update model grid parameters
model.Re     = Re;
model.long   = xq1 - ((max(xq1) + min(xq1))/2);
model.latg   = xq2 - ((max(xq2) + min(xq2))/2);
model.rzg    = (1000*Re) + xq3;
model.wlon   = max(model.long) - min(model.long);
model.wlat   = max(model.latg) - min(model.latg);
model.minrz  = min(model.rzg);
end