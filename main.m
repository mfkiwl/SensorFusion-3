%% Main runner
% Run this file from dataprocessing to optimization

clear; close all; clc;
%% Load Raw Data
imu = load('imu.mat');
gnss = load('gnss.mat');
lane = load('lane.mat');
can = load('can.mat');
snap = load('snap_raw.mat');

%% Pre Process raw data
dataset = dataprocessor(imu,gnss,can,lane,snap);
dataset.process();
% dataset.visualize();

%% Optimizer Options 

imu_ = dataset.proc_data.imu;
gnss_ = dataset.proc_data.gnss;
lane_ = dataset.proc_data.lane;
can_ = dataset.proc_data.can;
bias_ = dataset.proc_data.initBias;
t_ = dataset.proc_data.t;

% Covariance Setttings: Prior, Accel & Gyro Noise, Bias randomwalk 
covs_ = struct();
% Prior Covariance (Anchoring)
covs_.prior = struct();
covs_.prior.R = diag([0.1^2, 0.1^2, 0.1^2]);
covs_.prior.V = diag([1^2, 1^2, 1^2]);
covs_.prior.P = diag([5^2, 5^2, 5^2]);
covs_.prior.Bg = diag([1e-4, 1e-4, 1e-4]);
covs_.prior.Ba = diag([0.1^2, 0.1, 0.1^2]);
covs_.prior.WSF = 1e-2;
% covs_.prior.Params = diag([1e-4,1e-6,1e-4]);
covs_.prior.Params = 1e-5;

% IMU Accel, Gyro noise and bias randomwalk, ScaleFactorNoise
covs_.imu = struct();
covs_.imu.GyroscopeBiasNoise = 5e-10 * eye(3);
covs_.imu.GyroscopeNoise = 1e-5 * eye(3);
covs_.imu.AccelerometerBiasNoise = 5e-7* eye(3);
covs_.imu.AccelerometerNoise = 5/3 * 1e-3 * eye(3);
covs_.imu.ScaleFactorNoise = 1e-2;

% WSS 
covs_.wss = 1e-6 * eye(3);

% Optimization Options
options = struct();
options.CostThres = 1e-6;
options.StepThres = 1e-6;
options.IterThres = 50;
options.Algorithm = 'GN';
% GN : Gauss-Newton (Recommended for fast convergence, may not be stable for severely non-linear cases)
% LM : Levenberg-Marquardt(Not recommended for batch-wise optimization: wrong convergence)
% TR : Trust-Region (Recommended for stable convergence, but typically much slower than G-N method)

% If selected algorithm is LM, need to define parameters additionally
options.LM = struct();
options.LM.eta = 0.1;
options.LM.Lu = 11; % Lambda Up multiplier
options.LM.Ld = 5; % Lambda Decrease divider

% If selected algorithm is TR, need to define parameters additionally
options.TR = struct();
options.TR.eta1 = 0.6;
options.TR.eta2 = 0.9;
options.TR.gamma1 = 0.1;
options.TR.gamma2 = 2;
options.TR.thres = 1e-6; % Trust Region Radius Threshold

lane_.prev_num = 6; % Set preview number
lane_.prob_thres = 0.6; % Set lane prob threshold for discarding low-reliability data

%% INS + GNSS Fusion 
% sol = struct();
% sol.basic = optimizer(imu_,gnss_,lane_,can_,snap,bias_,t_,covs_,'basic',options);
% sol.basic.optimize();
% sol.basic.visualize();

%% INS + GNSS + WSS Fusion
% sol = struct();
% sol.partial = optimizer(imu_,gnss_,lane_,can_,snap,bias_,t_,covs_,'partial',options);
% sol.partial.optimize();
% sol.partial.visualize();

%% INS + GNSS + WSS + Lane Fusion

% To find warning id of most recent warning, use 'warning('query','last')'
warning('off','MATLAB:nearlySingularMatrix');

sol = struct();
% Optimize with INS + GNSS + WSS Fusion first 
sol.full = optimizer(imu_,gnss_,lane_,can_,snap,bias_,t_,covs_,'partial',options);
sol.full.optimize();
% Switch optimization mode to 2-phase and optimize with lane data
sol.full.opt.options.Algorithm = 'TR';
sol.full.update('2-phase'); % Update mode to 2-phase
sol.full.optimize();

%%

sol.full.visualize();

%% 
sol.full.map.visualize2DMap();


%%
sol.full.optimize2();