%% Parameter File for SynTomo
% All parameters required to run SynTomo are defined in this file. This
% file is the only input for SynTomo. See SynTomo manual for further
% description of parameters.
%
% B. VanderBeek (OCT-2020)
%

% Directories
msource = '../../toolbox'; % Matlab source-code for tomography program
taupjar = '../../TauP-2.4.5/lib/TauP-2.4.5.jar'; % Location of TauP jar-file
tomoin  = '../TomoInput/TEST'; % Input directory for tomography program
tomout  = '../TomoOutput/TEST'; % Output directory for tomography program

% General Parameters
% Forward Problem
npool            = 4; % Number of workers in Matlab parallel pool (npool = 0 will run serially)
ModelType        = 'P'; % Model type, 'P' or 'S'
RefModel         = 'iasp91'; % A reference 1D radial Earth model (see TauP manual for options)
tf_aniso         = true; % If true, includes anisotropy in travel-time calculation
% Inverse Problem
tf_invert        = true; % If true, perform inversion for isotroic velocity. If false, only predicts travel-times
tf_event_statics = true; % If true, include event static parameters in inverse problem. Generally, should be true.
tf_save_system   = false; % If true, will save system of equations defining tomography problem
Ldamp            = 1; % Damping multiplier for constrained least-squares inversion
Lsmooth          = [100,100,100]; % Smoothing multipliers for constrained least-squares inversion
lsqr_its         = 1e6; % Maximum number of allowed LSQR iterations (see lsqr doc).

% Forward Model Parameters
theModel       = 'model_SinkingBlock'; % Filename for existing or to-be-created model structure (stored in tomoin)
tf_build_model = true; % If true, will build model structure using below parameters and save to [tomoin,'/',theModel]
theVIZTOMO     = '../drex/SinkingBlock_Spherical/syntomo0020.h5'; % Location of VIZTOMO model
tf_cart        = false; % If true, the VIZTOMO model is defined in cartesian coordinates. If false, spherical coordinates assumed.
N1             = 223; % Number of nodes in longitudinal-direction
N2             = 223; % Number of nodes in latitudinal-direction
N3             = 67; % Number of nodes in radial-direction
order          = [1,3,2]; % Order of dimensions in VIZTOMO model

% Station Parameters
theStation       = 'station_SinkingBlock'; % Filename for existing or to-be-created station structure (stored in tomoin)
tf_build_station = true; % If true, will build station structure using below parameters and save to [tomoin,'/',theStation]
S1               = 25; % Number of stations in longitudinal-direction
S2               = 25; % Number of stations in latitudinal-direction
S1min            = -8.1479; % Minimum station longitude (degrees)
S1max            = 8.1479; % Maximum station longitude (degrees)
S2min            = -8.1479; % Minimum station latitude (degrees)
S2max            = 8.1479; % Maximum station latitude (degrees)

% Event Parameters
theEvent       = 'event_SinkingBlock'; % Filename for existing or to-be-created event structure (stored in tomoin)
tf_build_event = true; % If true, will build event structure using below parameters and save to [tomoin,'/',theEvent]
EAZIM          = linspace(0,330,12); % Vector of event azimuths with respect to model origin (deg. clockwise of north)
EDELT          = linspace(50,80,2); % Vector of event ranges (deg.)
EDPTH          = linspace(50,50,1); % Vector of event depths (km positive into Earth)

% Arrival Parameters
theArrival       = 'arrival_SinkingBlock_Ray_Iso'; % Filename for existing or to-be-created arrival structure (stored in tomoin)
tf_build_arrival = true; % If true, will build arrival structure using below parameters and save to [tomoin,'/',theArrival]
T                = 0; % Vector of observation periods (s; use 0 for ray theory modelling)
ERR              = 0.1; % A constant error associated with each arrival
aPhase           = 'P'; % Either 'P' or 'S'. Can only auto-generate arrival file for first arriving teleseismic phases

% Perturbational/Inverse Model Parameters
thePert       = 'pert_SinkingBlock'; % Filename for existing or to-be-created model perturbational structure (stored in tomoin)
tf_build_pert = true; % If true, will build perturbational structure using below parameters and save to [tomoin,'/',thePert]
M1            = 112; % Number of nodes in longitudinal-direction
M2            = 112; % Number of nodes in latitudinal-direction
M3            = 34; % Number of nodes in radial-direction
