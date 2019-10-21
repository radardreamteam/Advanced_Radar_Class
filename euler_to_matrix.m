% FUNCTIONNAME Short description of function
%
% Purpose
%
% inputs:
%
% outputs:
%

% References:
% [1]
%
% Notes:
%
% Author:
function R=euler_to_matrix(yaw_psi,pitch_theta,roll_phi)
%function R=euler_to_matrix(psi,the,phi)
% inputs:
% psi - yaw angle (rad):    rotation about current z-axis
% the - pitch angle (rad):  rotation about current y-axis
% phi - roll angle (rad):   rotation about current x-axis
%
% outputs:
% R - rotation matrix
%

% Mike Grabbe
% July 2004
N = length(yaw_psi);
if N>1
    
else
    r_psi=[cos(yaw_psi) -sin(yaw_psi) 0;sin(yaw_psi) cos(yaw_psi) 0;0 0 1];
    r_the=[cos(pitch_theta) 0 sin(pitch_theta);0 1 0;-sin(pitch_theta) 0 cos(pitch_theta)];
    r_phi=[1 0 0;0 cos(roll_phi) -sin(roll_phi);0 sin(roll_phi) cos(roll_phi)];
    
    R=r_psi*r_the*r_phi;
end