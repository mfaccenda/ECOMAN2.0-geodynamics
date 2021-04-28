function [tt,rayP,inc,toff,outPhase,iwarn] = taup_time(aModel,aPhase,delta,sdepth,rdepth)
% TAUP_PATH: A very simple matlab wrapper to call the taup_time function
% from the TauP Java toolbox.
%
% INPUT
%     aModel: A TauP recognized model name, e.g. 'iasp91' (string)
%     aPhase: A TauP recognized seismic phase name (string)
%      delta: A distance (deg.)
%     sdepth: Source depth (positive into Earth; km)
%     rdepth: Receiver depth (positive into Earth; km)
%
% OUTPUT
%         tt: Travel-time (s)
%       rayP: Ray parameter (s/deg)
%        inc: Ray incidence angle at station (deg. from vertical)
%       toff: Ray take-off angle at source (deg. from vertical)
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

% Create the travel-time java object. Takes ~1 ms to make.
ttobj = javaObject('edu.sc.seis.TauP.TauP_Time',aModel);

% Set the TauP_Time parameters
ttobj.setPhaseNames({aPhase});
ttobj.setSourceDepth(sdepth);
ttobj.setReceiverDepth(rdepth);

% Calculate arrival times
ttobj.calculate(delta);

% Extract results
N = ttobj.getNumArrivals;

% Check if arrival was returned
if (N == 0) && (strcmpi(aPhase,'P') || strcmpi(aPhase,'S'))
    % Find first arrival that matches input phase
    warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg. Looking for first ',aPhase,' arrival...']);
    if strcmpi(aPhase,'P')
        ttobj.setPhaseNames(phase_set(1,:));
    elseif strcmpi(aPhase,'S')
        ttobj.setPhaseNames(phase_set(2,:))
    end
    
    % Calculate arrival times
    ttobj.calculate(delta);
    
    % Select only the first arrival (if any returned)
    N = ttobj.getNumArrivals;
    N = min(N,1);
end

iwarn = 0;
if N > 0
    if N > 1
        iwarn = 2;
        warning(['Multiple ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.; using first.']);
    end
    ttout    = ttobj.getArrival(0); % Indexing starts at 0
    outPhase = char(ttout.getName);
    tt       = ttout.getTime;
    rayP     = ttout.getRayParam*pi/180;
    inc      = ttout.getIncidentAngle;
    toff     = ttout.getTakeoffAngle;
else
    iwarn = 1;
    warning(['No ',aPhase,' arrivals predicted at delta = ',num2str(delta),' deg.']);
    tt       = Inf;
    rayP     = Inf;
    inc      = Inf;
    toff     = Inf;
    outPhase = aPhase;
end
