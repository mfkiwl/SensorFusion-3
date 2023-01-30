%% Main runner
% Run this file from dataprocessing to optimization

clear; close all; clc;
%% Load Raw Data (Not 'really' raw, take a look at 'dataloader.py')

% [Sejong Dataset Scenarios]
% Scenario 1: In the current folder : Challenging Scenario
% Scenario 2: 2022-08-05--02-41-03 : plot with HD map (BRT)
% Scenario 3: 2022-08-05--02-57-24 : plot with HD map (세종청사)
% Scenario 4: 2022-08-05--03-13-16 : plot with HD map (세종청사)
% Scenario 5: 2022-08-05--03-45-16 : Convergence Error: almost same as sc4
% Scenario 6: 2022-08-05--04-19-33 : Challenging Scenario 2 (굳이?)

% base_path = 'D:\SJ_Dataset\2022-08-05\';
% scenario2 = '2022-08-05--02-41-03';
% scenario3 = '2022-08-05--02-57-24';
% scenario4 = '2022-08-05--03-13-16';

imu1 = load('imu.mat');
gnss1 = load('gnss.mat');
lane1 = load('lane.mat');
can1 = load('can.mat');
snap1 = load('snap_raw.mat');

% imu2 = load(strcat(base_path,scenario2,'\results\imu.mat'));
% gnss2 = load(strcat(base_path,scenario2,'\results\gnss.mat'));
% lane2 = load(strcat(base_path,scenario2,'\results\lane.mat'));
% can2 = load(strcat(base_path,scenario2,'\results\can.mat'));
% snap2 = load(strcat(base_path,scenario2,'\results\snap_raw.mat'));
% 
% imu3 = load(strcat(base_path,scenario3,'\results\imu.mat'));
% gnss3 = load(strcat(base_path,scenario3,'\results\gnss.mat'));
% lane3 = load(strcat(base_path,scenario3,'\results\lane.mat'));
% can3 = load(strcat(base_path,scenario3,'\results\can.mat'));
% snap3 = load(strcat(base_path,scenario3,'\results\snap_raw.mat'));
% 
% imu4 = load(strcat(base_path,scenario4,'\results\imu.mat'));
% gnss4 = load(strcat(base_path,scenario4,'\results\gnss.mat'));
% lane4 = load(strcat(base_path,scenario4,'\results\lane.mat'));
% can4 = load(strcat(base_path,scenario4,'\results\can.mat'));
% snap4 = load(strcat(base_path,scenario4,'\results\snap_raw.mat'));

% imu = load('imu.mat');
% gnss = load('gnss.mat');
% lane = load('lane.mat');
% can = load('can.mat');
% snap = load('snap_raw.mat');

% gnss_pos = [gnss1.pos(:,1),gnss1.pos(:,2);
%             gnss2.pos(:,1),gnss2.pos(:,2);
%             gnss3.pos(:,1),gnss3.pos(:,2);
%             gnss4.pos(:,1),gnss4.pos(:,2)];
% 
% snap_pos = [snap1.lat', snap1.lon';
%             snap2.lat', snap2.lon';
%             snap3.lat', snap3.lon';
%             snap4.lat', snap4.lon'];

figure(25);
% geoplot(snap_pos(:,1),snap_pos(:,2),'r.');
p1 = geoplot(snap1.lat,snap1.lon,'r--','LineWidth',1.5); hold on; grid on;
% p2 = geoplot(snap2.lat,snap2.lon,'g--','LineWidth',1.5);
% p3 = geoplot(snap3.lat,snap3.lon,'b--','LineWidth',1.5);
% % p4 = geoplot(snap4.lat,snap4.lon,'c--','LineWidth',1.5);

geobasemap satellite
title('Full Vehicle Trajectory for Sejong Dataset')
% legend([p1,p2,p3],'Scenario 1','Scenario 2','Scenario 3')

%% Pre Process raw data
imu = imu1; gnss = gnss1; lane = lane1; can = can1;
dataset = dataprocessor(imu,gnss,can,lane);
dataset.process();
% dataset.visualize();

imu_ = dataset.proc_data.imu;
gnss_ = dataset.proc_data.gnss;
lane_ = dataset.proc_data.lane;
can_ = dataset.proc_data.can;
bias_ = dataset.proc_data.initBias;
t_ = dataset.proc_data.t;
snap = [];

%% Optimizer Options 

% Covariance Setttings: Prior, Accel & Gyro Noise, Bias randomwalk 
covs_ = struct();
% Prior Covariance (Anchoring)
covs_.prior = struct();
covs_.prior.R = diag([0.1^2, 0.1^2, 0.1^2]);
covs_.prior.V = diag([1^2, 1^2, 1^2]);
covs_.prior.P = diag([5^2, 5^2, 5^2]);
covs_.prior.Bg = diag([1e-3, 1e-3, 1e-3]);
covs_.prior.Ba = diag([0.1^2, 0.1^2, 0.1^2]);
covs_.prior.WSF = 1e-4;
% covs_.prior.Params = diag([1e-4,1e-6,1e-4]);
covs_.prior.Params = 1e-5;

% IMU Accel, Gyro noise and bias randomwalk, ScaleFactorNoise
covs_.imu = struct();
covs_.imu.GyroscopeBiasNoise = 5e-9 * eye(3);
covs_.imu.GyroscopeNoise = 1e-5 * eye(3);
covs_.imu.AccelerometerBiasNoise = 5e-7* eye(3);
covs_.imu.AccelerometerNoise = 5/3 * 1e-3 * eye(3);
covs_.imu.ScaleFactorNoise = 1e-4;

% WSS 
% covs_.wss = diag([1e-2,1e-5,1e-5]);
covs_.wss = 3e-2 * eye(3);

% Optimization Options
options = struct();
options.CostThres = 1e-6;
options.StepThres = 1e-6;
options.IterThres = 500;
options.Algorithm = 'TR';
% GN : Gauss-Newton (Recommended for fast convergence, may not be stable for severely non-linear cases)
% TR : Trust-Region (Recommended for stable convergence, but typically much slower than G-N method)

% If selected algorithm is TR, need to define parameters additionally
options.TR = struct();
options.TR.eta1 = 0.5;
options.TR.eta2 = 0.9;
options.TR.gamma1 = 0.1;  
options.TR.gamma2 = 2;
options.TR.thres = 1e-6; % Trust Region Radius Threshold

lane_.prev_num = 6; % Set preview number
lane_.prob_thres = 0.8; % Set lane prob threshold for discarding low-reliability data
lane_.std_thres = 0.1; % Set lane std threshold for discarding unstable data
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
%%
sol.full.optimize();
% Switch optimization mode to 2-phase and optimize with lane data
sol.full.opt.options.Algorithm = 'TR';
%%
sol.full.update('2-phase'); % Update mode to 2-phase


%%
sol.full.map.dummyF();
% initArcParams = sol.full.map.arc_segments; % save data

%%
sol.full.map.dummy2();
%%
sol.full.map.dummy3();


%%

% To test phase 2 optimization, use saved parameters from phase 1
% sol.full.map.arc_segments = load('initArcParams.mat');

sol.full.optimize();

%%
% sol.full.map.validate();


sol.full.plotConfEllipse(0.9);

%%
sol.full.visualize();

%% 
sol.full.map.visualize2DMap();
% saved_map = sol.full.map; %% Save data 
%% Compare with HD Map
% HD Map Data obtained from: http://map.ngii.go.kr/ms/pblictn/preciseRoadMap.do
%
% cd D:/SJ_Dataset/HDMap/Map1/HDMap_UTMK_타원체고/
% cd D:/SJ_Dataset/HDMap/Map2/SEC01_BRT_내부간선도로/HDMap_UTMK_타원체고/
% cd D:/SJ_Dataset/HDMap/Map2/SEC02_세종정부청사_주변/HDMap_UTMK_타원체고/


% D:\SJ_Dataset\HDMap\Map1\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp
% D:\SJ_Dataset\HDMap\Map2\SEC01_BRT_내부간선도로\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp
% D:\SJ_Dataset\HDMap\Map2\SEC02_세종정부청사_주변\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp

% Load Saved data

% ArcMap1 = load("ResultData\MapSC1.mat");
% ArcMap2 = load("ResultData\MapSC2.mat");
ArcMap3 = load("ResultData\MapSC3.mat");
% ArcMap1.saved_map.visualize2DMap();
% ArcMap2.saved_map.visualize2DMap();
% ArcMap3.saved_map.visualize2DMap();

T1 = readgeotable("D:\SJ_Dataset\HDMap\Map2\SEC01_BRT_내부간선도로\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp");
T2 = readgeotable("D:\SJ_Dataset\HDMap\Map2\SEC02_세종정부청사_주변\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp");
% mapshow(T1); hold on;
% mapshow(T2);

sol.full.visualizeHD(T1,T2);
% sol.full.visualizeHD(T1,T2,ArcMap3.saved_map);

% for i=1:9
%     string = num2str(i);
%     front = "D:\SJ_Dataset\HDMap\Map3\";
%     back = "\HDMap_UTMK_타원체고\B2_SURFACELINEMARK.shp";
%     T1 = readgeotable(strcat(front,string,back));
%     mapshow(T1); hold on; grid on; axis equal;
% end
% back2 = "\HDMap_UTMK_타원체고\A1_NODE.shp";
% T2 = readgeotable(strcat(front,string,back));

% T1 = readgeotable("D:\NaverLabsDataset\pangyo\road_layout\1.0.0\pangyo_A3_LINK_3D.shp");
% mapshow(T1)

%% Test
T3 = readgeotable("D:\SJ_Dataset\HDMap\Map2\SEC01_BRT_내부간선도로\HDMap_UTMK_타원체고\A1_NODE.shp");
T4 = readgeotable("D:\SJ_Dataset\HDMap\Map2\SEC01_BRT_내부간선도로\HDMap_UTMK_타원체고\A2_LINK.shp");

%% 3D Plot for lanes
N = length(sol.full.states);
ZeroPreviewLeft = zeros(3,N);
ZeroPreviewRight = zeros(3,N);

for i=1:N
    R = sol.full.states{i}.R; P = sol.full.states{i}.P;
    laneL = sol.full.states{i}.left(:,1);
    laneR = sol.full.states{i}.right(:,1);
    ZeroPreviewLeft(:,i) = P + R * laneL;
    ZeroPreviewRight(:,i) = P + R * laneR;
end

figure(99);
plot3(ZeroPreviewLeft(1,:),ZeroPreviewLeft(2,:),ZeroPreviewLeft(3,:),'r.'); grid on; axis equal; hold on;
plot3(ZeroPreviewRight(1,:),ZeroPreviewRight(2,:),ZeroPreviewRight(3,:),'b.');

