function K = st_EvalHFFKernel(ray,model,pert,PhaseType,T,refModel,tf_aniso,XG,YG,ZG)
% ST_EVALHFFKERNEL: Computes the approximate finite-frequency travel-time
% through the model using the heuristic finite-frequency travel-time kernel
% (an approximation to the Born sensitivity kernel). Created for SynTomo
% matlab package.
%
% INPUT
%         ray: SynTomo raypath structure (see function 'st_GetRayPath1D')
%       model: SynTomo model structure (see SynTomo manual)
%        pert: SynTomo perturbational structure (see SynTomo manual)
%   PhaseType: Wave type. Use 'P' for compressional and 'S' for shear.
%           T: Dominant period of wave (s)
%    refModel: Name of reference 1D radial Earth model (e.g. 'iasp91'; see
%              TauP manual for recognized model names).
%    tf_aniso: If true, include anisotropy in travel-time computation.
%          XG: Global cartesian x-coordinates of model grid
%          YG: Global cartesian y-coordinates of model grid
%          ZG: Global cartesian z-coordinates of model grid
%
% OUTPUT
%   K: Kernel structure (see 'st_GetHFFKernel'). This function creates the 
%      following fields in 'K'...
%      ttff: Approximate finite-frequency travel-time through input model (s)
%      ttrt: Ray-theoretical travel-time through input model (s)
%      tt1D: TauP predicted travel-time through 1D model 'refModel'
%         w: Kernel travel-time sensitivity mapped to perturbational model
%            grid
%       ind: Linear index of weights in perturbational grid
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

% Ray Theoretical Travel-Times
if strcmp(PhaseType,'P')
    % Get average 1D velocities at ray path segments
    v1D  = get_TauP_model(abs(ray.rz),refModel);
    v1D  = (v1D(1:end-1) + v1D(2:end))./2;
    % Get anisotropic model velocities along ray path
    if tf_aniso
        % Average anisotropic velocities at ray segments
        vray = st_GetAnisoVel(model,ray.lon,ray.lat,ray.rz,ray.azm,ray.elv);
    else
        % Average isotropic velocities along ray segments
        vray = TrilinInterp(model.long,model.latg,model.rzg,model.Vp,ray.lon,ray.lat,ray.rz);
    end
    % Get 1D velocity profile through model
%     vrz = get_TauP_model(abs(model.zg),refModel);
    fld  = 'Vp';
elseif strcmp(PhaseType,'S')
    % Get average 1D velocities at ray path segments
    [~,v1D] = get_TauP_model(abs(ray.rz),refModel);
    v1D     = (v1D(1:end-1) + v1D(2:end))./2;
    % Average isotropic velocities along ray segments
    vray    = TrilinInterp(model.long,model.latg,model.rzg,ray.lon,ray.lat,ray.rz);
    % Get 1D velocity profile through model
%     [~,vrz] = get_TauP_model(abs(model.zg),refModel);
    fld      = 'Vs';
else
    error(['Unrecognized phase type ',PhaseType,'.']);
end

% Subset points for consideration in HFFK definition to reduce computation time 
[XG,YG,ZG,IND] = st_ReduceTestPoints(T,ray,vray,-model.wlon/2,model.wlon/2,...
    -model.wlat/2,model.wlat/2,model.Re,XG,YG,ZG);
% Get heuristic finite frequency kernel
K = st_GetHFFKernel(T,ray,vray,XG,YG,ZG,IND,mean(abs(diff(model.rzg)))/2);
% Map kernel indices to array subscripts
[ix,jy,kz] = ind2sub([model.n1,model.n2,model.n3],K.ind);

% Get apparent velocities at kernel nodes
if tf_aniso && strcmp(fld,'Vp')
    v3D = st_GetAnisoVel(model,model.long(ix),model.latg(jy),model.rzg(kz),K.azm,K.elv);
else
    v3D = model.(fld)(K.ind);
end

% Construct a 1D profile from the 3D model by averaging
% velocities in approximate kernel at each depth
nrz  = accumarray(kz,ones(size(K.ind)),[model.n3,1],[],1);
vrz  = accumarray(kz,model.(fld)(K.ind),[model.n3,1]);
vrz  = vrz./nrz;
% Interpolate 1D model to ray path segments
vref = interp1(model.rzg,vrz,ray.rz);
vref = (vref(1:end-1) + vref(2:end))./2;
vray = (vray(1:end-1) + vray(2:end))./2;
% Compute 1D and approximate finite frequency travel-times
dL     = diff(ray.d);
ttref  = sum(dL./vref);
K.ttff = ttref + sum(((1./v3D)-(1./vrz(kz))).*K.w);
K.ttrt = sum(dL./vray);
K.tt1D = sum(dL./v1D);

% Map kernel to perturbational grid
% Interpolation weights
[ix,jy,kz] = ind2sub([model.n1,model.n2,model.n3],K.ind);
[~,w,ind]  = TrilinInterp(pert.u.long,pert.u.latg,pert.u.rzg,[],model.long(ix),model.latg(jy),model.rzg(kz));
% Accumulate kernel values on perturbational nodes
K.w    = w.*repmat(K.w,1,8);
K.w    = K.w(:);
ind    = ind(:);
K.ind  = unique(ind);
[~,iw] = ismember(ind,K.ind);
K.w    = accumarray(iw,K.w,[max(iw),1]);
