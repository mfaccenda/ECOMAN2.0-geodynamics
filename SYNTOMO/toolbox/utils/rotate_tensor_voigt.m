function Cr = rotate_tensor_voigt(C,alpha,beta,gamma,order)
% ROTATE_TENSOR_VOIGT: Rotates a 6x6 tensor in Voigt notation.
%
% INPUT
%       C: A 6x6 tensor in Voigt notation.
%   alpha: Rotation about x1-axis; positive clockwise looking at origin
%          from positive end of x1-axis (degrees).
%    beta: Rotation about x2-axis; positive clockwise looking at origin
%          from positive end of x2-axis (degrees).
%   gamma: Rotation about x3-axis; positive clockwise looking at origin 
%          from positive end of x3-axis (degrees).
%   order: A 3-element vector defining the order in which to apply the
%          rotations (e.g. 'order = [3,2,1]' would rotate first around the
%          x3-axis by gamma, then the x2-axis by beta, lastly around the
%          x1-axis by alpha).
%
% OUTPUT
%      Cr: The rotated tensor.
%
% NOTES
% + This routine is a simplified version of the function 'MS_rot3' found in
%   the Matlab Seismic Anisotropy (MSAT) Toolbox (Walker & Wookey, 2012).
%
% REFERENCES
% Walker, A. M., & Wookey, J. (2012). MSAT—A new toolkit for the analysis 
%   of elastic and seismic anisotropy. Computers & Geosciences, 49, 81-90.
%
% B. VanderBeek (OCT-2020)
%

% Radian angles
alpha = alpha*pi/180;
beta  = beta*pi/180;
gamma = gamma*pi/180;

% Construct 3x3 rotation matrix
R        = zeros(3,3,3) ;
R(1,:,:) = [1, 0, 0; 0, cos(alpha), sin(alpha); 0, -sin(alpha), cos(alpha)];
R(2,:,:) = [cos(beta), 0, -sin(beta); 0, 1, 0; sin(beta), 0, cos(beta)];
R(3,:,:) = [cos(gamma), sin(gamma), 0; -sin(gamma), cos(gamma), 0; 0, 0, 1];
R        = squeeze(R(order(3),:,:))*squeeze(R(order(2),:,:))...
           *squeeze(R(order(1),:,:));

% Form the K matrix (based on Bowers 'Applied Mechanics of Solids', Chapter 3)
K1 = [R(1,1).^2, R(1,2).^2, R(1,3).^2; ...
      R(2,1).^2, R(2,2).^2, R(2,3).^2; ...
      R(3,1).^2, R(3,2).^2, R(3,3).^2];

K2 = [R(1,2).*R(1,3), R(1,3).*R(1,1), R(1,1).*R(1,2); ...
      R(2,2).*R(2,3), R(2,3).*R(2,1), R(2,1).*R(2,2); ...
      R(3,2).*R(3,3), R(3,3).*R(3,1), R(3,1).*R(3,2)];

K3 = [R(2,1).*R(3,1), R(2,2).*R(3,2), R(2,3).*R(3,3); ...
      R(3,1).*R(1,1), R(3,2).*R(1,2), R(3,3).*R(1,3); ...
      R(1,1).*R(2,1), R(1,2).*R(2,2), R(1,3).*R(2,3)] ;

K4 = [R(2,2).*R(3,3)+R(2,3).*R(3,2),...
      R(2,3).*R(3,1)+R(2,1).*R(3,3),...
      R(2,1).*R(3,2)+R(2,2).*R(3,1); ...
      R(3,2).*R(1,3)+R(3,3).*R(1,2),...
      R(3,3).*R(1,1)+R(3,1).*R(1,3),...
      R(3,1).*R(1,2)+R(3,2).*R(1,1); ...
      R(1,2).*R(2,3)+R(1,3).*R(2,2),...
      R(1,3).*R(2,1)+R(1,1).*R(2,3),...
      R(1,1).*R(2,2)+R(1,2).*R(2,1)];

K = [K1, 2*K2; ...
     K3,   K4];

% Apply rotation
Cr = K*C*(K');
