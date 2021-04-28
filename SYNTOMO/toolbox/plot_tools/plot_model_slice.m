function [C,H] = plot_model_slice(xg,yg,zg,F,slice,varargin)
% PLOT_MODEL_SLICE: Simple function to plot cross-section through
% seismic model. Created for SynTomo matlab package.
%
% INPUT
%      xg: x1-coordinate vector of model array
%      yg: x2-coordinate vector of model array
%      zg: x3-coordinate vector of model array
%       F: Model array
%   slice: For a horizontal cross-section plane, 'slice' specifies depth of
%          cross-section. For a vertically oriented cross-section plane,
%          'slice' is a 1x3 vector defining the x,y-origin and azimuth
%          (degrees CCW of x-axis) of the cross-section plane.
% <optional>
% + These are name-value-pair arguments specified as
%   plot_model_slice(...,'name',value). Options for 'name' are...
%    'tf_cont': Logical that if true will draw a countour plot. Default is
%               false.
%       'cint': Array of contour intervals if 'tf_cont' is true.
%       'cmap': A Nx3 colormap array (e.g. cmap = jet)
%       'limc': A 1x2 array defining the minimum and maximum of colorscale
%       'limx': A 1x2 array defining plot limits in x1-direction
%       'limy': A 1x2 array defining plot limits in x2-direction
%
% OUTPUT
%   C: Cross-section structure (see get_cross_section)
%   H: Handle to figure
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

% Arrange plot options for easy reference
popts = reshape(varargin,2,length(varargin)/2)';

% Get the cross-section plane
C = get_cross_section(xg,yg,zg,F,slice);
bx = [C.X1(:,1);C.X1(end,:)';flipud(C.X1(:,end));fliplr(C.X1(1,:))'];
by = [C.X2(:,1);C.X2(end,:)';flipud(C.X2(:,end));fliplr(C.X2(1,:))'];

% Identify plot type
[~,iopt] = ismember('tf_cont',popts(:,1));
if iopt > 0
    tf_cont = popts{iopt,2};
else
    tf_cont = false;
end

% Make plot
H = figure;
hold on;
box on;
if tf_cont
    % Check for contour intervals
    [~,iopt] = ismember('cint',popts(:,1));
    if iopt > 0
        G = contourf(C.X1,C.X2,C.F,popts{iopt,2});
    else
        G = contourf(C.X1,C.X2,C.F);
    end
else
    G = pcolor(C.X1,C.X2,C.F);
    set(G,'EdgeColor','none');
end
% Plot outline of cross-section
plot(bx,by,'-k','linewidth',1);

% Assign colormap
[~,iopt] = ismember('cmap',popts(:,1));
if iopt > 0
    colormap(popts{iopt,2});
end

% Assign axes limits
[~,iopt] = ismember('limx',popts(:,1));
if iopt > 0
    xlim(popts{iops,2});
end
[~,iopt] = ismember('limy',popts(:,1));
if iopt > 0
    xlim(popts{iops,2});
end

% Assign color limits
[~,iopt] = ismember('limc',popts(:,1));
if iopt > 0
    caxis(popts{iopt,2});
end

% Derive title
if length(slice) == 3
    tag = ['x_i = ',num2str(slice(1)),'; y_i = ',num2str(slice(2)),'; azimuth = ',num2str(slice(3))];
else
    tag = ['z_i = ',num2str(slice)];
end
title(tag);
colorbar;
axis image;

