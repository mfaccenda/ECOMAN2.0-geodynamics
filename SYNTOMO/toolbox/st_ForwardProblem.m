function [misfit,Jac] = st_ForwardProblem(event,station,arrival,model,pert,...
    npool,refModel,PhaseType,tf_aniso,tf_event_static)
% ST_FORWARDPROBLEM: SynTomo function that controls the computation of
% travel-times and partial derivatives for teleseismic arrivals. Created 
% for SynTomo Matlab package.
%
% INPUT
%             event: SynTomo event structure defining teleseismic event
%                    parameters (see SynTomo manual).
%           station: SynTomo station structure defining seismic array
%                    (see SynTomo manual).
%           arrival: SynTomo arrival structure defining travel-time
%                    observations (see SynTomo manual).
%             model: SynTomo model structure defining seismic velocity
%                    and/or elastic model used in forward modeling arrivals
%                    (see SynTomo manual).
%              pert: SynTomo perturbational model structure defining the inversion
%                    inversion grid (see SynTomo manual).
%             npool: Number of workers in Matlab pool for parallel
%                    computations (use 0 for serial computations/debugging)
%          refModel: Name of reference 1D radial Earth model (e.g.
%                    'iasp91'; see TauP manual for recognized model names).
%         PhaseType: Wave type; use 'P' for compressional and 'S' for shear
%          tf_aniso: If true, includes anisotropy in the forward modelling
%                    of travel-times. If false, uses isotropic model.
%   tf_event_static: If true, adds event static terms to the Jacobian
%                    matrix. This is a common practice in regional 
%                    teleseismic tomography when structure outside the
%                    array is unknown.
%
% OUTPUT
%   misfit: SynTomo structure containing information regarding the forward
%           calculation with the following fields...
%           runtime_forward: Elapsed real time spent in forward calculation (s)
%                   ttime1D: TauP predicted travel-time through 'refModel' (s)
%                     ttime: Travel-time predicted through input model (s)
%                       res: Travel-time residual (i.e. arrival.time - misfit.ttime; s)
%                      lray: Length of ray path in model (km)
%                       rms: The root-mean-square of the travel-time residuals (s)
%                      wrms: The error-weighted root-mean-square of the travel-time residuals (s)
%                      chi2: The chi-squared value of the residuals (--)
%      Jac: Sparse Jacobian matrix filled with travel-time derivatives.
%           Rows correspond to a travel-time observation and columns
%           correspond to a perturbational model parameter.
%
% NOTES
% + All travel-times are absolute source-to-receiver times
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

% Allow common indexing of event and station coordinates for parallel loop
[~,ievt] = ismember(arrival.event,event.id);
[~,ista] = ismember(arrival.station,station.name);
elat     = event.latitude(ievt);
elon     = event.longitude(ievt);
edpt     = event.depth(ievt);
slat     = station.latitude(ista);
slon     = station.longitude(ista);


% Define sliced structures for parallel loop
narr = length(arrival.ttime);
% Arrival structure
arrvSlice(narr).event   = [];
arrvSlice(narr).station = [];
arrvSlice(narr).phase   = [];
arrvSlice(narr).period  = [];
arrvSlice(narr).ttime   = [];
for ii = 1:narr
    arrvSlice(ii).event   = arrival.event(ii);
    arrvSlice(ii).station = arrival.station(ii);
    arrvSlice(ii).phase   = arrival.phase{ii};
    arrvSlice(ii).period  = arrival.period(ii);
    arrvSlice(ii).ttime   = arrival.ttime(ii);
    arrvSlice(ii).error   = arrival.error(ii);
end

% Jacobian structure
Jac(narr).row = [];
Jac(narr).col = [];
Jac(narr).val = [];

% Predicted data and data residual (misfit) structure
misfitSlice(narr).ttime1D = [];
misfitSlice(narr).ttime   = [];
misfitSlice(narr).res     = [];
misfitSlice(narr).lray    = [];

% Some more variable parsing to make the parfor loop happy
nparam    = pert.u.n1*pert.u.n2*pert.u.n3;
Re        = model.Re;
minlon    = -model.wlon/2;
maxlon    = model.wlon/2;
minlat    = -model.wlat/2;
maxlat    = model.wlat/2;
minz      = model.minrz;
maxz      = 0;
dz        = (maxz - minz)/(model.n3 - 1);

% Global cartesian grid for finite-frequency calculations. We only want to
% make these once since this mapping can take a non-trivial amount of time.
[YG,XG,ZG] = meshgrid(model.latg,model.long,model.rzg);
[XG,YG,ZG] = geo2cart(XG,YG,Re + ZG,0,0,Re,0);

% Forward Problem--Compute travel-times and partial derivatives
fprintf(['\n Solving forward problem for ',num2str(narr),' arrivals.\n']);
t1 = tic;
parfor (ii = 1:narr, npool)
    % Get 1D ray path within model domain via TauP
    ray = st_GetRayPath1D(refModel,arrvSlice(ii).phase,elat(ii),elon(ii),edpt(ii),...
        slat(ii),slon(ii),[minlat,maxlat],[minlon,maxlon],[minz,maxz],dz,Re);
    
    % Recompute travel-times using true model velocities
    if arrvSlice(ii).period == 0
        K    = st_EvalRayKernel(ray,model,pert,PhaseType,refModel,tf_aniso);
        tt1D = K.tt1D;
        tt3D = K.tt3D;
    else
        K    = st_EvalHFFKernel(ray,model,pert,PhaseType,arrvSlice(ii).period,refModel,tf_aniso,XG,YG,ZG);
        tt1D = K.tt1D;
        tt3D = K.ttff;
    end
    
    % Fill misfit structure
    misfitSlice(ii).ttime1D = ray.tt_taup;
    misfitSlice(ii).ttime   = ray.tt_taup + (tt3D - tt1D);
    misfitSlice(ii).res     = arrvSlice(ii).ttime - misfitSlice(ii).ttime;
    misfitSlice(ii).lray    = ray.d(end);
    
    %% Fill Jacobian
    
    % Include event statics?
    if tf_event_static
        K.ind = cat(1,K.ind,nparam + ievt(ii));
        K.w   = cat(1,K.w,1);
    end
    
    % Store Jacobian values
    Jac(ii).row = ii*ones(length(K.ind),1);
    Jac(ii).col = K.ind;
    Jac(ii).val = K.w./arrvSlice(ii).error;
end
t2 = toc(t1);
fprintf(['\n Time spent in forward problem: ',num2str(t2),' s.\n']);

% Store computation time
misfit.runtime_forward = t2;

% Concatonate misfit structure
misfit.ttime1D = cat(1,misfitSlice(:).ttime1D);
misfit.ttime   = cat(1,misfitSlice(:).ttime);
misfit.res     = cat(1,misfitSlice(:).res);
misfit.lray    = cat(1,misfitSlice(:).lray);

% Model fit metrics
misfit.rmsr  = rms(misfit.res);
misfit.wrmsr = sqrt(sum((misfit.res./arrival.error).^2)/sum(1./(arrival.error.^2)));
misfit.chi2  = mean((misfit.res./arrival.error).^2);

% Construct sparse Jacobian
iJ   = cat(1,Jac(:).row);
jJ   = cat(1,Jac(:).col);
vJ   = cat(1,Jac(:).val);
Jac  = sparse(iJ,jJ,vJ);
