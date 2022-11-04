%%  Test ArcFit Jacobian
clear; close all; clc;

% 
n = 10; % number of sub-segments
syms x0 y0 tau0 xd yd
kappa = sym('kappa',[1,n]);
L = sym('L',[1,n]);

%% Propagate


heading = tau0;

% for i=1:n
%     if i == 1
%         xc = x0 - 1/kappa(i) * sin(heading);
%         yc = y0 + 1/kappa(i) * cos(heading);
%     else
%         xc = xc + (1/kappa(i-1) - 1/kappa(i)) * sin(heading);
%         yc = yc - (1/kappa(i-1) - 1/kappa(i)) * cos(heading);
%     end
%     heading = heading + kappa(i) * L(i);
% end
% 
% res = sqrt((xc - xd)^2 + (yc - yd)^2) - abs(1/kappa(n));
% 
% param_jac = jacobian(res,[x0,y0,tau0]);
% kappa_jac = jacobian(res,kappa);
% L_jac = jacobian(res,L);


disp(jacobian(abs(1/kappa(1)),kappa(1)))