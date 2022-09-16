clear; close all; clc;
%% Load Raw Data
imu = load('imu.mat');
gnss = load('gnss.mat');
lane = load('lane.mat');
can = load('can.mat');

%% Pre Process raw data
dataset = dataprocessor(imu,gnss,can,lane);
dataset.process();
% dataset.visualize();

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
covs_.prior.Params = diag([0.1^2,1e-6,0.1^2]);

% IMU Accel, Gyro noise and bias randomwalk, ScaleFactorNoise
covs_.imu = struct();
covs_.imu.GyroscopeBiasNoise = 5e-10 * eye(3);
covs_.imu.GyroscopeNoise = 1e-5 * eye(3);
covs_.imu.AccelerometerBiasNoise = 5e-6* eye(3);
covs_.imu.AccelerometerNoise = 5/3 * 1e-3 * eye(3);
covs_.imu.ScaleFactorNoise = 1e-6;

% WSS 
covs_.wss = 1e-4 * eye(3);

% Optimization Options
options = struct();
options.CostThres = 1e-6;
options.StepThres = 1e-6;
options.IterThres = 30;
options.Algorithm = 'GN';
% GN : Gauss-Newton
% LM : Levenberg-Marquardt(To Be Done)
% TR : Trust-Region(To Be Done)



%% INS + GNSS Fusion 
% sol = struct();
sol.basic = optimizer(imu_,gnss_,lane_,can_,bias_,t_,covs_,'basic',options);
sol.basic.optimize();
%%
sol.basic.visualize();

%% INS + GNSS + WSS Fusion
% sol = struct();
% sol.partial = optimizer(imu_,gnss_,lane_,can_,bias_,t_,covs_,'partial',options);
% sol.partial.optimize();
% sol.partial.visualize();
