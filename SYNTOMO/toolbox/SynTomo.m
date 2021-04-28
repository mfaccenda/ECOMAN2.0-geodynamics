function SynTomo(theParFile)
%% SYNTOMO: Main function that runs the SynTomo package.
%
% INPUT
%   theParFile: Name of Matlab parameter file for running SynTomo (see
%               SynTomo manual for description).
%
% OUTPUT
%   Results of forward and inverse modelling are stored in the directories
%   specified in 'theParFile'.
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

% Set parameter file
run(theParFile);

% Force the TauP path to be a full path. I think the relative paths require
% a call to 'javaaddpath' every time TauP is called which substantially
% increases the runtime.
taupjar = dir(taupjar); %#ok The 'taupjar' variable is created above by 'run(theParFile)'.
taupjar = [taupjar.folder,'/',taupjar.name];

% Add paths
addpath(genpath(msource),'-begin');
setenv('TAUPJAR',taupjar);
% Check that input/output directories exist. If not, create them.
if ~isfolder(tomoin)
    system(['mkdir -p ',tomoin]);
end
if ~isfolder(tomout)
    system(['mkdir -p ',tomout]);
end

% Start parallel pool
pp = gcp('nocreate');
if isempty(pp)
    if npool > 0
        parpool('local',npool,'AttachedFiles',{getenv('TAUPJAR')});
    end
else
    fprintf('\n Found existing parallel pool. Re-initializing parallel pool.\n');
    if (npool > 0) && ~(pp.NumWorkers == npool)
        delete(gcp);
        parpool('local',npool,'AttachedFiles',{getenv('TAUPJAR')});
    else
        delete(gcp);
    end
end

% Track total run time
ttot = tic;

%% Build/Load Inputs

% Save a copy of the parameter file in the ouput directory for reference
system(['cp ',theParFile,' ',tomout,'/ParamSynTomo.m']);

% Model
if tf_build_model
    fprintf('\n Building model structure...\n');
    model = st_ReadModel(theVIZTOMO,N1,N2,N3,tf_cart,order);
    save([tomoin,'/',theModel,'.mat'],'model');
else
    load([tomoin,'/',theModel,'.mat'],'model');
end
% If only doing isotropic modelling, remove elastic tensor fields to
% reduce model size
if ~tf_aniso && isfield(model,'RHO')
    model = rmfield(model,{'RHO','C11','C12','C13','C14','C15','C16',...
        'C22','C23','C24','C25','C26','C33','C34','C35','C36','C44',...
        'C45','C46','C55','C56','C66'});
end

% Station
if tf_build_station
    fprintf('\n Building station structure...\n');
    % Define evenly spaced stations
    kmpdeg      = model.Re*pi/180;
    Lx          = kmpdeg*(S1max - S1min);
    Ly          = kmpdeg*(S2max - S2min);
    Xsta        = linspace(-Lx/2,Lx/2,S1);
    Ysta        = linspace(-Ly/2,Ly/2,S2);
    % Create regular array of stations
    [Ysta,Xsta] = meshgrid(Ysta,Xsta);
    Xsta        = Xsta(:);
    Ysta        = Ysta(:);
    % Map stations back to geographic coordinates
    D           = sqrt((Xsta.^2) + (Ysta.^2))./kmpdeg;
    AZ          = atan2d(Xsta,Ysta);
    AZ          = wrapTo360(AZ);
    [Ysta,Xsta] = reckon(0,0,D,AZ);
    % Create station structure (define station names and store coordinates)
    station.name = cell(length(Xsta),1);
    fmt = ['%0',num2str(1 + floor(log10(length(Xsta)))),'i'];
    for ii = 1:length(Xsta)
        station.name{ii} = ['SYN',num2str(ii,fmt)];
    end
    station.longitude = Xsta;
    station.latitude  = Ysta;
    save([tomoin,'/',theStation,'.mat'],'station');
else
    load([tomoin,'/',theStation,'.mat'],'station');
end

% Event
if tf_build_event
    fprintf('\n Building event structure...\n');
    % Create array of events
    [edel,eazim,edpth] = meshgrid(EDELT,EAZIM,EDPTH);
    eazim = eazim(:);
    edel  = edel(:);
    edpth = edpth(:);
    % Create event structure
    event.id    = (1:length(edpth))';
    event.depth = edpth(:);
    [event.latitude,event.longitude] = reckon(0,0,edel,eazim);
    % Save
    save([tomoin,'/',theEvent,'.mat'],'event');
else
    load([tomoin,'/',theEvent,'.mat'],'event');
end

% Arrival
if tf_build_arrival
    fprintf('\n Building dummy arrival structure...\n');
    nsta            = length(station.name);
    nevt            = length(event.id);
    arrival.event   = reshape(repmat(event.id(:)',nsta,1),nsta*nevt,1);
    arrival.station = reshape(repmat(station.name(:),1,nevt),nsta*nevt,1);
    arrival.phase   = repmat({aPhase},nsta*nevt,1);
    arrival.period  = T(1)*ones(nsta*nevt,1);
    iT = 2;
    while iT <= length(T)
        arrival.event   = cat(1,arrival.event,arrival.event);
        arrival.station = cat(1,arrival.station,arrival.station);
        arrival.phase   = cat(1,arrival.phase,arrival.phase);
        arrival.period  = cat(1,arrival.period,T(iT)*ones(nsta*nevt,1));
    end
    arrival.ttime   = zeros(length(arrival.period),1);
    arrival.error   = ERR*ones(length(arrival.period),1);
    % Do not need to save dummy arrival file
else
    load([tomoin,'/',theArrival,'.mat'],'arrival');
end

% Perturbational Model
if tf_build_pert
    fprintf('\n Building perturbational model structure...\n');
    pert.u.long = linspace(-model.wlon/2,model.wlon/2,M1)';
    pert.u.latg = linspace(-model.wlat/2,model.wlat/2,M2)';
    pert.u.rzg  = linspace(0,model.minrz,M3)';
    pert.u.n1   = length(pert.u.long);
    pert.u.n2   = length(pert.u.latg);
    pert.u.n3   = length(pert.u.rzg);
    % Save
    save([tomoin,'/',thePert,'.mat'],'pert');
else
    load([tomoin,'/',thePert,'.mat'],'pert');
end

%% Call Tomography Routines

% Forward problem
[misfit,Jac] = st_ForwardProblem(event,station,arrival,model,pert,...
    npool,RefModel,ModelType,tf_aniso,tf_event_statics);

% Define synthetic arrival structure and 1D starting model and then
% recompute travel-times for inversion using 1D model.
if ~tf_invert
    % If not doing an inversion, save an arrival structure with the
    % predicted times
    fprintf('\n Creating synthetic arrival structure...\n');
    arrival.ttime = misfit.ttime;
    % Save
    save([tomoin,'/',theArrival,'_SYN.mat'],'arrival');
elseif tf_build_arrival && tf_invert
    % If we have to build the arrival structure and invert, then we must
    % also create a 1D model
    fprintf('\n Creating synthetic arrival structure...\n');
    arrival.ttime = misfit.ttime;
    fprintf('\n Creating 1D model structure and re-computing travel-times for inversion...\n');
    % Make model 1D
    fldnm = fieldnames(model);
    for ifld = 1:length(fldnm)
        [n1,n2,n3] = size(model.(fldnm{ifld}));
        if (n1 == model.n1) && (n2 == model.n2) && (n3 == model.n3)
            model.(fldnm{ifld}) = repmat(model.(fldnm{ifld})(1,1,:),model.n1,model.n2,1);
        end
    end
    % Save
    save([tomoin,'/',theArrival,'_SYN.mat'],'arrival');
    save([tomoin,'/',theModel,'_i1j1START1D.mat'],'model');
    
    % Re-run forward problem
    [misfit,Jac] = st_ForwardProblem(event,station,arrival,model,pert,...
        npool,RefModel,ModelType,tf_aniso,tf_event_statics);
end

% Stop the parallel pool
if npool > 0
    delete(gcp);
end

% Perform inversion?
if tf_invert
    % Inverse problem
    [pert,misfit,event,model,D,Sx,Sy,Sz] = st_InverseProblem(Jac,misfit,...
        arrival,pert,model,event,tf_event_statics,ModelType,Ldamp,Lsmooth,lsqr_its);
    % Save results
    save([tomout,'/','pert.mat'],'pert');
    save([tomout,'/','misfit.mat'],'misfit');
    save([tomout,'/','model.mat'],'model');
    if tf_event_statics
        save([tomout,'/','event.mat'],'event');
    end
    % Save Jacobian and constraint equations
    if tf_save_system
        save([tomout,'/','Jac.mat'],'Jac');
        save([tomout,'/','D.mat'],'D');
        save([tomout,'/','Sx.mat'],'Sx');
        save([tomout,'/','Sy.mat'],'Sy');
        save([tomout,'/','Sz.mat'],'Sz');
    end
else
    % Save predicted travel-times
    save([tomout,'/','misfit.mat'],'misfit');
    % Save Jacobian
    if tf_save_system
        save([tomout,'/','Jac.mat'],'Jac');
    end
end
ttot = toc(ttot);
fprintf(['\n Total runtime for SynTomo: ',num2str(ttot),' s.\n']);