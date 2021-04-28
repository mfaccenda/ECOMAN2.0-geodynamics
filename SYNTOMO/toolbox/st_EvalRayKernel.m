function K = st_EvalRayKernel(ray,model,pert,PhaseType,refModel,tf_aniso)
% ST_EVALRAYKERNEL: Computes the ray-theoretical travel-time through the 
% input model and sensitivity to perturbational model parameters. Created 
% for SynTomo matlab package.
%
% INPUT
%         ray: SynTomo raypath structure (see function 'st_GetRayPath1D')
%       model: SynTomo model structure (see SynTomo manual)
%        pert: SynTomo perturbational structure (see SynTomo manual)
%   PhaseType: Wave type. Use 'P' for compressional and 'S' for shear.
%    refModel: Name of reference 1D radial Earth model (e.g. 'iasp91'; see
%              TauP manual for recognized model names).
%    tf_aniso: If true, include anisotropy in travel-time computation.
%
% OUTPUT
%   K: Ray theoretical sensitivity kernel structure. This function creates 
%      the following fields in 'K'...
%      tt1D: TauP predicted travel-time through 1D model 'refModel' (s)
%      tt3D: Ray-theoretical travel-time through input model (s)
%         w: Ray-theoretical travel-time sensitivity kernel mapped to the
%            perturbational model grid
%       ind: Linear index of weights in perturbational grid
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
        vray = (vray(1:end-1) + vray(2:end))./2;
    else
        % Average isotropic velocities along ray segments
        vray = TrilinInterp(model.long,model.latg,model.rzg,model.Vp,ray.lon,ray.lat,ray.rz);
        vray = (vray(1:end-1) + vray(2:end))./2;
    end
elseif strcmp(PhaseType,'S')
    % Get average 1D velocities at ray path segments
    [~,v1D] = get_TauP_model(abs(ray.rz),refModel);
    v1D     = (v1D(1:end-1) + v1D(2:end))./2;
    % Average isotropic velocities along ray segments
    vray    = TrilinInterp(model.long,model.latg,model.rzg,ray.lon,ray.lat,ray.rz);
    vray    = (vray(1:end-1) + vray(2:end))./2;
else
    error(['Unrecognized phase type ',PhaseType,'.']);
end
% Ray theoretical travel-time
dL     = diff(ray.d);
K.tt1D = sum(dL./v1D);
K.tt3D = sum(dL./vray);

% Map slowness partials to perturbational grid
% Get linear interpolation weights
[~,w,ind] = TrilinInterp(pert.u.long,pert.u.latg,pert.u.rzg,[],ray.lon,ray.lat,ray.rz);
% Scale weights by ray segment lengths
dL = abs(diff(ray.d));
dL = cat(1,dL(1),(dL(1:end-1) + dL(2:end)),dL(end))./2;
w  = w.*repmat(dL,1,8);
% Accumulate weights
ind    = ind(:);
w      = w(:);
K.ind  = unique(ind);
[~,iw] = ismember(ind,K.ind);
K.w    = accumarray(iw,w,[max(iw),1]);
