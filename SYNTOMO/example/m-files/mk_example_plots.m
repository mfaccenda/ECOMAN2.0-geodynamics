%% Make Example Plots
% This script gives a few examples of how to plot results from SynTomo.
close all;
clear;
clc;

% SynTomo Inputs
theStation  = '../TomoInput/TEST/station_SinkingBlock.mat';
theEvent    = '../TomoInput/TEST/event_SinkingBlock.mat';
theArrival  = '../TomoInput/TEST/arrival_SinkingBlock_Ray_Iso_SYN.mat';
theStart    = '../TomoInput/TEST/model_SinkingBlock_i1j1START1D.mat';
theTrue     = '../TomoInput/TEST/model_SinkingBlock.mat';
% SynTomo Outputs
theMisfit   = '../TomoOutput/TEST/misfit.mat';
theModel    = '../TomoOutput/TEST/model.mat';
thePert     = '../TomoOutput/TEST/pert.mat';

%% Plot Mean Station Delays
% Provides coarse view of vertically averaged structure beneath stations

% Load data
load(theStation,'station');
load(theEvent,'event');
load(theArrival,'arrival');
load(theMisfit,'misfit');
load(theStart,'model');

% Associate arrivals to event and stations
[~,ievt] = ismember(arrival.event,event.id);
[~,ista] = ismember(arrival.station,station.name);

% Demean data by event
DT  = misfit.res;
MED = accumarray(ievt,DT,[length(event.id),1],@mean);
DT  = DT - MED(ievt);

% Mean station delays
MSD = accumarray(ista,DT,[length(station.name),1],@mean);

% Plot
figure;
scatter(station.longitude,station.latitude,100,MSD,'filled');
colormap(jet); colorbar; caxis([-1 1]);
axis image; box on; grid on;
xlim([-model.wlon,model.wlon]./2); ylim([-model.wlat,model.wlat]./2);
title('Mean Station Delays');
xlabel('longitude (deg.)'); ylabel('latitude (deg.)');

%% Plot Cross-sections Through Model

% Cross-section parameters
x0  = 0; % x-origin of cross-section
y0  = 0; % y-origin of cross-section
az  = 0; % azimuth of cross-section (deg. CCW of east)
zx  = -300; % Depth for cross-section in XY-plane (km)
fld = 'Vp'; % Field in model to plot
% Plot parameters
tf_cont = false; % If true, plot contours
cint    = -0.1:0.005:0.1; % Contour levels if above is true
cmap    = flipud(jet); % Colormap
limc    = [-0.025 0.025]; % Colormap limits
tf_mask = false; % If true, mask result based on derivative weight sum
dws_min = 10; % Mask model where derivative weight sum is less than dws_min

% Define fractional perturbations
% True model
load(theStart,'model');
V0 = model.(fld); clear('model');
load(theTrue,'model');
Ftrue = (model.(fld) - V0)./V0;
% Recovered model
clear('model');
load(theModel,'model');
F = (model.(fld) - V0)./V0;
% Apply mask
load(thePert,'pert');
DWS = interp3(pert.u.latg(:)',pert.u.long(:),pert.u.rzg(:)',pert.u.DWS,model.latg(:)',model.long(:),model.rzg(:)');
if tf_mask
    F(DWS < dws_min) = NaN;
end

% Plot data coverage (vertical cross-section)
plot_model_slice(pert.u.long,pert.u.latg,pert.u.rzg,pert.u.DWS,[x0,y0,az],...
    'tf_cont',false,'cmap',hot(64),'limc',[0,500]);
title('DWS'); xlabel('km'); ylabel('km');

% Plot data coverage (depth cross-section)
plot_model_slice(pert.u.long,pert.u.latg,pert.u.rzg,pert.u.DWS,zx,...
    'tf_cont',false,'cmap',hot(64),'limc',[0,500]);
title('DWS'); xlabel('km'); ylabel('km');

% Plot true perturbations (vertical cross-section)
plot_model_slice(model.long,model.latg,model.rzg,Ftrue,[x0,y0,az],...
    'tf_cont',tf_cont,'cmap',cmap,'limc',limc,'cint',cint);
title('True dlnV'); xlabel('km'); ylabel('km');

% Plot recovered perturbations (vertical cross-section)
plot_model_slice(model.long,model.latg,model.rzg,F,[x0,y0,az],...
    'tf_cont',tf_cont,'cmap',cmap,'limc',limc,'cint',cint);
title('Recovered dlnV'); xlabel('km'); ylabel('km');

% Plot true perturbations (depth section)
plot_model_slice(model.long,model.latg,model.rzg,Ftrue,zx,...
    'tf_cont',tf_cont,'cmap',cmap,'limc',limc,'cint',cint);
title('True dlnV'); xlabel('km'); ylabel('km');

% Plot recovered perturbations (depth section)
[C,H] = plot_model_slice(model.long,model.latg,model.rzg,F,zx,...
    'tf_cont',tf_cont,'cmap',cmap,'limc',limc,'cint',cint);
title('Recovered dlnV'); xlabel('km'); ylabel('km');
