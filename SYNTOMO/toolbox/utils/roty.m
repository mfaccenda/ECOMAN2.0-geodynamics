function R = roty(beta)
% ROTX: Returns rotation matrix for a rotation about y-axis. Rotation angle
% 'beta' is measured counterclockwise positive looking at origin from
% positive y-axis in degrees.
%
% B. VanderBeek (OCT-2020)
%
% Created to reproduce results of function of same name in Matlab's Phase
% Array System toolbox which is not always included in the Matlab
% installation.
%

R = [cosd(beta), 0, sind(beta);...
     0,          1,          0;...
    -sind(beta), 0, cosd(beta)];