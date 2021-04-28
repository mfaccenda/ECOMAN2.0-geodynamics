function [pert,misfit,event,model,D,Sx,Sy,Sz] = st_InverseProblem(Jac,...
    misfit,arrival,pert,model,event,tf_event_static,PhaseType,Ldamp,Lsmooth,lsqr_its)
% ST_INVERSEPROBLEM: Function that defines and solves the constrained 
% least-squares tomography problem. Created for SynTomo Matlab package.
%
% INPUT
%               Jac: Jacobian matrix composed of travel-time derivatives
%                    defined by 'st_ForwardProblem'.
%            misfit: SynTomo structure containing travel-time residual
%                    information (see SynTomo manual).
%           arrival: SynTomo arrival structure defining travel-time
%                    observations (see SynTomo manual).
%              pert: SynTomo perturbational model structure defining the
%                    inversion grid (see SynTomo manual).
%             model: SynTomo model structure defining seismic velocity
%                    and/or elastic model used in forward modeling arrivals
%                    (see SynTomo manual).
%             event: SynTomo event structure defining teleseismic event
%                    parameters (see SynTomo manual)
%   tf_event_static: If true, event statics are included as inversion
%                    parameters.
%         PhaseType: Specifies compressional ('P') or shear ('S') wave 
%                    speed inversion.
%             Ldamp: Scalar multiplier for inversion damping constraint
%           Lsmooth: A 1x3 vector defining scalar multipliers for the
%                    smoothing constraints in the x-, y-, and z-directions.
%          lsqr_its: Number of iterations allowed for the LSQR algorithm
%
% OUTPUT
%     pert: SynTomo perturbational model structure updated with the
%           following fields...
%             vprior: Starting model velocities interpolated to the 
%                     perturbational model.
%                 dm: Slowness perturbations
%                DWS: Derivative weight sum
%            dm_norm: The norm of the slowness perturbations
%           dm_rough: The roughness of the slowness perturbations (i.e.
%                     norm of the vector defined by the multiplication of 
%                     the perturbational vector with the array defining the
%                     Laplacian smoothing coefficients.
%   misfit: SynTomo misfit structure updated with the following fields...
%           runtime_inverse: Total time required to perform inversion (s)
%                  lsqr_its: Number of LSQR iterations performed
%                   res_est: Estimate of travel-time residual vector with
%                            respect to updated velocity model
%                  rmsr_est: Root-mean-square of estimated travel-time
%                            residuals
%                 wrmsr_est: Weighted root-mean-square of estimated
%                            travel-time residuals
%                  chi2_est: Chi-squared value of estimated travel-time
%                            residuals
%    event: SynTomo event structure with the field 'static' added that
%           contains the event static parameters determined by the
%           inversion
%    model: SynTomo model structure with updated wave speeds
%        D: Sparse array of damping equations
%       Sx: Sparse array of Laplacian smoothing coefficients in x-direction
%       Sy: Sparse array of Laplacian smoothing coefficients in y-direction
%       Sz: Sparse array of Laplacian smoothing coefficients in z-direction
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

% Define some parameters
nparam = pert.u.n1*pert.u.n2*pert.u.n3;
if tf_event_static
    nstat = length(event.id);
else
    nstat = 0;
end

if strcmp(PhaseType,'P')
    vfld = 'Vp';
elseif strcmp(PhaseType,'S')
    vfld = 'Vs';
else
    error(['Unrecognized phase type ',PhaseType])
end

% Initial model velocities at perturbational nodes
[Yq,Xq,Zq] = meshgrid(pert.u.latg,pert.u.long,pert.u.rzg);
Vi = TrilinInterp(model.long,model.latg,model.rzg,model.(vfld),Xq(:),Yq(:),Zq(:));
pert.u.vprior = reshape(Vi,pert.u.n1,pert.u.n2,pert.u.n3);


% Construct Sparse Damping Constraint
fprintf('\n Building smoothing and damping equations.\n');
id = (1:nparam)';
wd = Ldamp*Vi;
D  = sparse(id,id,wd,nparam,nparam + nstat);


% Construct Laplacian smoothing constraints
[Sx,isx,jsx] = LaplaceCoeff(pert.u.long,1,[pert.u.n1,pert.u.n2,pert.u.n3]);
[Sy,isy,jsy] = LaplaceCoeff(pert.u.latg,2,[pert.u.n1,pert.u.n2,pert.u.n3]);
[Sz,isz,jsz] = LaplaceCoeff(pert.u.rzg,3,[pert.u.n1,pert.u.n2,pert.u.n3]);
% Scale Laplacian coefficients
Sx = Lsmooth(1)*Vi(jsx).*Sx;
Sy = Lsmooth(2)*Vi(jsy).*Sy;
Sz = Lsmooth(3)*Vi(jsz).*Sz;
% Build sparse smoothing matrices
Sx = sparse(isx,jsx,Sx,nparam,nparam + nstat);
Sy = sparse(isy,jsy,Sy,nparam,nparam + nstat);
Sz = sparse(isz,jsz,Sz,nparam,nparam + nstat);


% Solve inverse problem
fprintf('\n Solving inverse problem via LSQR.\n');
t1 = tic;
[dm,~,~,it] = lsqr(cat(1,Jac,D,Sx,Sy,Sz),...
    cat(1,misfit.res./arrival.error,zeros(4*nparam,1)),[],lsqr_its);
t2 = toc(t1);
fprintf(['\n Time spent in LSQR: ',num2str(t2),' s.\n']);
misfit.runtime_inverse = t2;
misfit.lsqr_its        = it;

% Parse solution
if tf_event_static
    event.static = dm((nparam+1):end);
end
pert.u.dm = reshape(dm(1:nparam),pert.u.n1,pert.u.n2,pert.u.n3);

% Store derivative weight sum
pert.u.DWS = full(sum(Jac,1));
pert.u.DWS = reshape(pert.u.DWS(1:nparam),pert.u.n1,pert.u.n2,pert.u.n3);

% Store solution norm and roughness
pert.u.dm_norm   = norm(dm(:));
pert.u.dm_rough  = norm((Sx./Lsmooth(1))*dm(:))...
    + norm((Sy./Lsmooth(2))*dm(:))...
    + norm((Sz./Lsmooth(3))*dm(:));

% Store estimate of solution residual/fit
misfit.res_est   = misfit.res - (Jac*dm(:)).*arrival.error;
misfit.rmsr_est  = rms(misfit.res_est);
misfit.wrmsr_est = sqrt(sum((misfit.res_est./arrival.error).^2)/sum(1./(arrival.error.^2)));
misfit.chi2_est  = mean((misfit.res_est./arrival.error).^2);

% Update the forward model
dm = interp3(pert.u.latg(:)',pert.u.long(:),pert.u.rzg(:)',pert.u.dm,...
    model.latg(:)',model.long(:),model.rzg(:)');
model.(vfld) = 1./((1./model.(vfld)) + dm);
