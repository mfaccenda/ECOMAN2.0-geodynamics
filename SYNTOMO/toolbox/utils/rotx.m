function R = rotx(alpha)
% ROTX: Returns rotation matrix for a rotation about x-axis. Rotation angle
% 'alpha' is measured counterclockwise positive looking at origin from
% positive x-axis in degrees.
%
% B. VanderBeek (OCT-2020)
%
% Created to reproduce results of function of same name in Matlab's Phase
% Array System toolbox which is not always included in the Matlab
% installation.
%

R = [1,           0,            0;...
     0, cosd(alpha), -sind(alpha);...
     0, sind(alpha), cosd(alpha)];