function R = rotz(gamma)
% ROTX: Returns rotation matrix for a rotation about z-axis. Rotation angle
% 'gamma' is measured counterclockwise positive looking at origin from
% positive z-axis in degrees.
%
% B. VanderBeek (OCT-2020)
%
% Created to reproduce results of function of same name in Matlab's Phase
% Array System toolbox which is not always included in the Matlab
% installation.
%

R = [cosd(gamma), -sind(gamma), 0;...
     sind(gamma),  cosd(gamma), 0;...
     0,           0,            1];