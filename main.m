%% Main runner
% Run this file from dataprocessing to optimization

clear; close all; clc;
%% Load Raw Data (Not 'really' raw, take a look at 'dataloader.py')

% [Sejong Dataset Scenarios]
% Scenario 1: In the current folder : Challenging Scenario
% Scenario 2: 2022-08-05--02-41-03 : plot with HD map (BRT)
% Scenario 3: 2022-08-05--02-57-24 : plot with HD map (세종청사)
% Scenario 4: 2022-08-05--03-13-16 : dataprocessing error --> need to fix
% Scenario 5: 2022-08-05--03-45-16 : dataprocessing error --> need to fix
% Scenario 6: 2022-08-05--04-19-33 : Challenging Scenario 2 (굳이?)

base_path = 'D:\SJ_Dataset\2022-08-05\';
scenario = '2022-08-05--04-19-33';

imu = load(strcat(base_path,scenario,'\results\imu.mat'));
gnss = load(strcat(base_path,scenario,'\results\gnss.mat'));
lane = load(strcat(base_path,scenario,'\results\lane.mat'));
can = load(strcat(base_path,scenario,'\results\can.mat'));

% imu = load('imu.mat');
% gnss = load('gnss.mat');
% lane = load('lane.mat');
% can = load('can.mat');
% snap = load('snap_raw.mat');

%% Pre Process raw data
dataset = dataprocessor(imu,gnss,can,lane);
dataset.process();
% dataset.visualize();

%% Optimizer Options 

imu_ = dataset.proc_data.imu;
gnss_ = dataset.proc_data.gnss;
lane_ = dataset.proc_data.lane;
can_ = dataset.proc_data.can;
bias_ = dataset.proc_data.initBias;
t_ = dataset.proc_data.t;
snap = 0;

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
covs_.imu.ScaleFactorNoise = 1e-6;

% WSS 
covs_.wss = 1e-4 * eye(3);

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
lane_.minL = 5; % Minimum arc length (to prevent singularity)

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
%%
sol.full.update('2-phase'); % Update mode to 2-phase


%%
sol.full.map.dummyF();
% initArcParams = sol.full.map.arc_segments; % save data

%%
sol.full.map.dummy.initFit.visualize();


%%

% To test phase 2 optimization, use saved parameters from phase 1
% sol.full.map.arc_segments = load('initArcParams.mat');

sol.full.optimize();

%%
% sol.full.map.validate();




%%
sol.full.visualize();

%% 
sol.full.map.visualize2DMap();

%% Compare with HD Map
% HD Map Data obtained from: http://map.ngii.go.kr/ms/pblictn/preciseRoadMap.do
%
% cd D:/SJ_Dataset/HDMap/Map1/HDMap_UTMK_타원체고/
% cd D:/SJ_Dataset/HDMap/Map2/SEC01_BRT_내부간선도로/HDMap_UTMK_타원체고/
% cd D:/SJ_Dataset/HDMap/Map2/SEC02_세종정부청사_주변/HDMap_UTMK_타원체고/


% D:\SJ_Dataset\HDMap\Map1\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp
% D:\SJ_Dataset\HDMap\Map2\SEC01_BRT_내부간선도로\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp
% D:\SJ_Dataset\HDMap\Map2\SEC02_세종정부청사_주변\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp

T1 = readgeotable("D:\SJ_Dataset\HDMap\Map1\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp");
sol.full.visualizeHD(T1);

