function [latray,lonray,zray,ttray,dray,outPhase] = taup_path(aModel,aPhase,slat,slon,sdepth,rlat,rlon)
% TAUP_PATH: A very simple matlab wrapper to call the taup_path function
% from the TauP Java toolbox. Returns a 1D ray path.
%
% INPUT
%     aModel: A TauP recognized model name, e.g. 'iasp91' (string)
%     aPhase: A TauP recognized seismic phase name (string)
%       slat: Source latitude (dec. deg.)
%       slon: Source longitude (dec. deg.)
%     sdepth: Source depth (positive into Earth; km)
%       rlat: Receiver latitude (dec. deg.)
%       rlon: Receiver longitude (dec. deg.)
%
% OUTPUT
%     latray: Latitude of ray path (dec. deg.)
%     lonray: Longitude of ray path (dec. deg.)
%       zray: Ray path depth (negative into the Earth with surface at 0; km)
%      ttray: Travel-time along ray path (s)
%       dray: Great circle distance along ray path (deg.)
%   outPhase: Output phase type. Should match aPhase. However, if aPhase is
%             'P' or 'S', then the first arriving compressional or shear
%             phase is returned.
%
% NOTES
% + The Java matlab object can compute a variety of ray path attributes
%   that may be added to this function if needed (e.g. pierce points).
%
% REFERENCES
% + Crotwell, H. P., T. J. Owens, and J. Ritsema (1999). The TauP ToolKit: 
%   Flexible Seismic Travel-Time and Raypath Utilities, Seismological 
%   Research Letters
%   http://www.seis.sc.edu
%
% B. VanderBeek AUG-2019
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


% Check if we need to add TauP jar. We don't want to always add it because
% it can be slow to add to the java class path.
jpath = javaclasspath('-dynamic');
if ~any(strcmp(jpath,getenv('TAUPJAR')))
    if isempty(getenv('TAUPJAR'))
        error('Missing ''TAUPJAR'' environment variable pointing to TauP jar file.');
    else
        javaaddpath(getenv('TAUPJAR'));
    end
end

% First arrival phase set
phase_set = {'P','Pdiff','PKIKP';'S','SKS','SKIKS'};

% Derive some values
[delta,az] = distance(slat,slon,rlat,rlon);

% Create the travel-time java object. Takes ~1 ms to make.
ttpath = javaObject('edu.sc.seis.TauP.TauP_Path',aModel);

% Set the TauP_Time parameters
ttpath.setPhaseNames({aPhase});
ttpath.setSourceDepth(sdepth);

% Calculate arrival times
ttpath.calculate(delta);

% Extract results
N = ttpath.getNumArrivals;

% Check if ray path was returned
if (N == 0) && (strcmp(aPhase,'P') || strcmp(aPhase,'S'))
    % Find first arrival that matches input phase
    warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg. Looking for first ',aPhase,' arrival...']);
    if strcmp(aPhase,'P')
        ttpath.setPhaseNames(phase_set(1,:));
    elseif strcmp(aPhase,'S')
        ttpath.setPhaseNames(phase_set(2,:))
    end
    
    % Calculate arrival times
    ttpath.calculate(delta);
    
    % Select only the first arrival (if any returned)
    N = ttpath.getNumArrivals;
    N = min(N,1);
end

% Define ray path
if N > 0
    if N > 1
        warning(['Multiple ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.; using first.']);
    end
    ttpath_out = ttpath.getArrival(0); % Indexing starts at 0
    outPhase   = ttpath_out.getName;
    Npath      = ttpath_out.getNumPathPoints;
    dray       = zeros(Npath,1);
    ttray      = zeros(Npath,1);
    zray       = zeros(Npath,1);
    for ii = 1:Npath
        P         = ttpath_out.getPathPoint(ii-1);
        dray(ii)  = P.getDistDeg;
        ttray(ii) = P.getTime;
        zray(ii)  = -P.getDepth;
    end
    % Get latitude/longitude coordinates of ray
    [latray,lonray] = reckon(slat,slon,dray,az);
else
    warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.']);
    latray   = [];
    lonray   = [];
    zray     = [];
    ttray    = [];
    dray     = [];
    outPhase = aPhase;
end
