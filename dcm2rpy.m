function rpy = dcm2rpy(R)
% dcm2rpy Converts Rotation Matrix to Roll, Pitch, Yaw angles 
% (Assumes Intrinsic Rotation whose Tait-Bryan angles are alpha, beta, gamma, about axes z-y-x) 
% 
% R = Rz(alpha) * Ry(beta) * Rx(gamma)
%   = [[ca * cb, ca * sb * sg - sa * cg, ca * sb * cg + sa * sg],
%      [sa * cb, sa * sb * sg + ca * cg, sa * sb * cg - ca * sg],
%      [-sb    , cb * sg               , cb * cg               ]]
% 
% alpha (yaw) = atan2(R(2,1),R(1,1))
% beta  (pitch) = -asin(R(3,1))
% gamma (roll) = atan2(R(3,2),R(3,3))
% 
% * Due to singularity for RPY or Euler angle representation, 
%   when plotting results, there may be jumping discontinuity

    phi = atan2(R(3,2),R(3,3));
    theta = -asin(R(3,1));
    psi = atan2(R(2,1),R(1,1));

    rpy = [phi;theta;psi];
end