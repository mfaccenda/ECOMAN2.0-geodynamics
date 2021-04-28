function [vp,vs,rho,Qp,Qs,Rz] = get_TauP_model(Rz,aModel)
% GET_TAUP_MODEL: Reads standard models included with the TauP Toolkit
%
% INPUT
%       Rz: Radial earth depth (km; positive into the Earth)
%           + If empty, returns original model file
%   aModel: A TauP model name (string)
%           + Assumes TauP model files are in Matlab search path
%           + Currently available models: iasp91, ak135, prem
%
% OUTPUT
%       vp: P-wave speed (km/s)
%       vs: S-wave speed (km/s)
%      rho: Density (g/cm^3)
%       Qp: Compressional quality factor
%       Qs: Shear quality factor
%
% NOTES
% + Easily add custom models by editing switch-case statement
% + This function may need to be updated if format of TauP standard models
%   changes.
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

% Find TauP models folder
pathtaup = getenv('TAUPJAR');
if ~isempty(pathtaup)
    pathtaup = strsplit(pathtaup,'/lib');
    pathtaup = pathtaup{1};
    pathtaup = [pathtaup,'/StdModels'];
else
    % Try to determine path to TauP models
    pathtaup = which('iasp91.tvel');
    if ~isempty(pathtaup)
        pathtaup = strsplit(pathtaup,'/iasp91');
        pathtaup = pathtaup{1};
    else
        error('Could not locate TauP standard models. Add to search path or define ''TAUPJAR'' environment variable pointing to TauP jar file.');
    end
end


switch aModel
    case {'iasp91','ak135'}
        % Read the file
        fid = fopen([pathtaup,'/',aModel,'.tvel']);
        M   = textscan(fid,'%f%f%f%f','Delimiter',' ','MultipleDelimsAsOne',true,'HeaderLines',2);
        fclose(fid);
        % Process
        M = cell2mat(M);
        n = size(M,1);
        M = cat(2,M,1e6*ones(n,2));
    case 'prem'
        % Read the file. Contains three named discontinuities
        fid = fopen([pathtaup,'/',aModel,'.nd']);
        M   = textscan(fid,'%f%f%f%f%f%f','Delimiter',' ','MultipleDelimsAsOne',true,'CommentStyle','m');
        M   = cat(1,M,textscan(fid,'%f%f%f%f%f%f','Delimiter',' ','MultipleDelimsAsOne',true,'CommentStyle','o'));
        M   = cat(1,M,textscan(fid,'%f%f%f%f%f%f','Delimiter',' ','MultipleDelimsAsOne',true,'CommentStyle','i'));
        fclose(fid);
        % Process
        M = cell2mat(M);
    otherwise
        error('Unrecognized model name.')
end

% Identify discontinuities and slightly perturb depths such that depth
% vector is monotonically increasing.
dz    = diff(M(:,1));
dz    = min(dz(dz > 0));
idisc = (diff(M(:,1)) == 0);
M([false;idisc],1) = M([false;idisc],1) + (dz/1000);

% Interpolate to desired model spacing
if ~isempty(Rz)
    vp  = interp1(M(:,1),M(:,2),Rz,'linear',M(end,2));
    vs  = interp1(M(:,1),M(:,3),Rz,'linear',M(end,3));
    rho = interp1(M(:,1),M(:,4),Rz,'linear',M(end,4));
    Qp  = interp1(M(:,1),M(:,5),Rz,'linear',M(end,5));
    Qs  = interp1(M(:,1),M(:,6),Rz,'linear',M(end,6));
else
    Rz  = M(:,1);
    vp  = M(:,2);
    vs  = M(:,3);
    rho = M(:,4);
    Qp  = M(:,5);
    Qs  = M(:,6);
end