clear; close all; clc;
%% Load data from previous research and convert format to usable one
load('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\VehicleLocalizationandMapping\output.mat');
load('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\VehicleLocalizationandMapping\lane.mat');

% Scenario 1 DSTART = 14800; DEND = 16000 --> 기준은 lane에 있는 t
% Scenario 2 DSTART = 17000; DEND = 21800

% Alter timestamp value for output.ubloxgps.t(79)
predT = output.ubloxgps.t(78) + 0.05;
delta = predT - output.ubloxgps.t(79);
output.ubloxgps.t(79:end) = output.ubloxgps.t(79:end) + delta;
% plot(output.ubloxgps.t(1:100) - output.ubloxgps.t(1));


%% Data Processing with DaeJeon Dataset
sc = 1;
if sc == 1
    DSTART = 14800; DEND = 16000;
elseif sc == 2
    DSTART = 17000; DEND = 21800;
end
% DSTART = 14800; DEND = 16000;
% DSTART = 17000; DEND = 21500;
gnss_raw = output.ubloxgps;
imu_raw = output.xs;
lane_raw = output.cam;
lla_ref = [mean(output.rtgps.lat),mean(output.rtgps.lon),mean(output.rtgps.alt)];

% Find GNSS Idx bnds
[~,gnssStartIdx] = min(abs(gnss_raw.t - lane.t(DSTART)));
[~,gnssEndIdx] = min(abs(gnss_raw.t - lane.t(DEND)));

sampledGNSSCov = gnss_raw.cov(:,gnssStartIdx:gnssEndIdx);
validIdxs = gnssStartIdx - 1 + find(sampledGNSSCov(1,:) < 2^2); % 2m threshold for GNSS outage filtering

% Find GNSS Outage intervals
for i=1:length(validIdxs)-1
    if validIdxs(i+1) - validIdxs(i) > 1
        disp(['GNSS Outage from Idx ',num2str(validIdxs(i)+1),' ~ ',num2str(validIdxs(i+1)-1)])
    end
end

state_t = [];
gnss_ = struct();
gnss_.state_idxs = [];
gnss_.hAcc = [];
gnss_.vAcc = [];
gnss_.pos = [];
gnss_.t = [];
gnss_.lla0 = lla_ref;
yaw = [];
for i=1:length(validIdxs)
    if i == length(validIdxs)
        state_t = [state_t, gnss_raw.t(validIdxs(i))];
        gnss_.t = [gnss_.t, gnss_raw.t(validIdxs(i))];
        
        [~,yawIdx] = min(abs(output.ubloxori.t - gnss_raw.t(validIdxs(i))));
        yaw = [yaw, output.ubloxori.rz(yawIdx)];

        gnss_.state_idxs = [gnss_.state_idxs, length(state_t)];
        gnss_.pos = [gnss_.pos; gnss_raw.lat(validIdxs(i)), gnss_raw.lon(validIdxs(i)), gnss_raw.lat(validIdxs(i))];
        gnss_.hAcc = [gnss_.hAcc, 9 * gnss_raw.cov(1,validIdxs(i))];
        gnss_.vAcc = [gnss_.vAcc, 9 * gnss_raw.cov(9,validIdxs(i))];
        break;
    end

    if validIdxs(i+1) - validIdxs(i) > 1 % Outage
%         disp(length(state_t)); error('1')
        gnss_.state_idxs = [gnss_.state_idxs, length(state_t)+1];
        state_t = [state_t, gnss_raw.t(validIdxs(i)):0.1:gnss_raw.t(validIdxs(i+1))];
        gnss_.t = [gnss_.t, gnss_raw.t(validIdxs(i))];
        gnss_.pos = [gnss_.pos; gnss_raw.lat(validIdxs(i)), gnss_raw.lon(validIdxs(i)), gnss_raw.lat(validIdxs(i))];
        gnss_.hAcc = [gnss_.hAcc, 9 * gnss_raw.cov(1,validIdxs(i))];
        gnss_.vAcc = [gnss_.vAcc, 9 * gnss_raw.cov(9,validIdxs(i))];

        [~,yawIdx] = min(abs(output.ubloxori.t - gnss_raw.t(validIdxs(i))));
        yaw = [yaw, output.ubloxori.rz(yawIdx)];

        % Fill in outage with 10Hz sampling time
        if abs(state_t(end) - gnss_raw.t(validIdxs(i+1))) < 1e-6
            state_t = state_t(1:end-1);
        end
    else
        state_t = [state_t, gnss_raw.t(validIdxs(i))];
        gnss_.t = [gnss_.t, gnss_raw.t(validIdxs(i))];
        gnss_.state_idxs = [gnss_.state_idxs, length(state_t)];
        gnss_.pos = [gnss_.pos; gnss_raw.lat(validIdxs(i)), gnss_raw.lon(validIdxs(i)), gnss_raw.lat(validIdxs(i))];
        gnss_.hAcc = [gnss_.hAcc, 9 * gnss_raw.cov(1,validIdxs(i))];
        gnss_.vAcc = [gnss_.vAcc, 9 * gnss_raw.cov(9,validIdxs(i))];

        [~,yawIdx] = min(abs(output.ubloxori.t - gnss_raw.t(validIdxs(i))));
        yaw = [yaw, output.ubloxori.rz(yawIdx)];
    end    
end
gnss_.vNED = [lane.vel.y(DSTART), lane.vel.x(DSTART), -lane.vel.z(DSTART)];
gnss_.bearing = 90 - 180/pi * lane.eul.z(DSTART);

% Cluster IMU values
% accel, gyro: n * 3 format 
% t
bias_ = struct();
bias_.accel = zeros(1,3);
bias_.gyro = zeros(1,3);
imu_ = [];

interpedAccelX = interp1(imu_raw.t,imu_raw.ax,state_t);
interpedAccelY = interp1(imu_raw.t,imu_raw.ay,state_t);
interpedAccelZ = interp1(imu_raw.t,imu_raw.az,state_t);
interpedGyroX = interp1(imu_raw.t,imu_raw.wx,state_t);
interpedGyroY = interp1(imu_raw.t,imu_raw.wy,state_t);
interpedGyroZ = interp1(imu_raw.t,imu_raw.wz,state_t);
base_t = state_t(1);

for i=1:length(state_t)-1
    t_lb = state_t(i); t_ub = state_t(i+1);
    imu = struct();    

    % Sample IMU
    sampledT = imu_raw.t(imu_raw.t >= t_lb & imu_raw.t < t_ub);
    sampledAccelX = imu_raw.ax(imu_raw.t >= t_lb & imu_raw.t < t_ub);
    sampledAccelY = imu_raw.ay(imu_raw.t >= t_lb & imu_raw.t < t_ub);
    sampledAccelZ = imu_raw.az(imu_raw.t >= t_lb & imu_raw.t < t_ub);
    sampledGyroX = imu_raw.wx(imu_raw.t >= t_lb & imu_raw.t < t_ub);
    sampledGyroY = imu_raw.wy(imu_raw.t >= t_lb & imu_raw.t < t_ub);
    sampledGyroZ = imu_raw.wz(imu_raw.t >= t_lb & imu_raw.t < t_ub);

    AccelX = [interpedAccelX(i), sampledAccelX];
    AccelY = [interpedAccelY(i), sampledAccelY];
    AccelZ = [interpedAccelZ(i), sampledAccelZ];
    GyroX = [interpedGyroX(i), sampledGyroX];
    GyroY = [interpedGyroY(i), sampledGyroY];
    GyroZ = [interpedGyroZ(i), sampledGyroZ];


    imu.accel = [AccelX', AccelY', AccelZ'];
    imu.gyro = [GyroX', GyroY', GyroZ'];
    imu.t = [t_lb - base_t, sampledT - base_t, t_ub - base_t];
    imu_ = [imu_, {imu}];
end

for i=1:length(imu_)
    if size(imu_{i}.accel,1) <= 1
        disp(['Less than 2 acceleration & gyro readings, please interpolate for stability: IMU Idx',num2str(i)])
    end
end

% Match RT values (reference)
ref = struct();
ref_t = output.rtgps.t;
ref.pos = [];

can_ = struct();
can_.state_idxs = 1:1:length(state_t);
can_.whl_spd = zeros(1,length(state_t));
for i=1:length(state_t)
    [~,ref_idx] = min(abs(ref_t - state_t(i)));
    ref.pos = [ref.pos; output.rtgps.lat(ref_idx), output.rtgps.lon(ref_idx), output.rtgps.alt(ref_idx)];
    
    % Compute Longitudinal Velocity
    th = output.rtimu.rz(ref_idx);
    R2d = [cos(th), -sin(th);
           sin(th), cos(th)];
    v_b = R2d' * [output.rtvel.vx(ref_idx), output.rtvel.vy(ref_idx)]';
    can_.whl_spd(i) = v_b(1);
end

% Process lane info
lane_ = lane_raw;
lane_.state_idxs = zeros(1,length(state_t));

for i=1:length(state_t)
    [~,Idx] = min(abs(lane_.t - state_t(i)));
    lane_.state_idxs(i) = Idx;
end
lane_.ly = lane_.l_inter;
lane_.lz = zeros(size(lane_.l_inter,1),11);

lane_.ry = lane_.r_inter;
lane_.rz = zeros(size(lane_.r_inter,1),11);

snap = ref.pos(1,:);
% snap = [];
t_ = state_t - base_t;

%% Optimizer Options 

% Covariance Setttings: Prior, Accel & Gyro Noise, Bias randomwalk 
covs_ = struct();
% Prior Covariance (Anchoring)
covs_.prior = struct();
covs_.prior.R = diag([0.1^2, 0.1^2, 0.1^2]);
covs_.prior.V = diag([2^2, 2^2, 2^2]);
covs_.prior.P = diag([1e-5, 1e-5, 1e-5]);
covs_.prior.Bg = diag([1e-3, 1e-3, 1e-3]);
covs_.prior.Ba = diag([0.3^2, 0.3^2, 0.3^2]);
covs_.prior.WSF = 1e-3;
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
covs_.wss = 1e-5 * eye(3);

% Optimization Options
options = struct();
options.CostThres = 1e-6;
options.StepThres = 1e-6;
options.IterThres = 500;
options.Algorithm = 'TR';
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
options.TR.eta1 = 0.5;
options.TR.eta2 = 0.9;
options.TR.gamma1 = 0.1;  
options.TR.gamma2 = 2;
options.TR.thres = 1e-6; % Trust Region Radius Threshold

lane_.prev_num = 6; % Set preview number
lane_.prob_thres = 0.6; % Set lane prob threshold for discarding low-reliability data
lane_.minL = 5; % Minimum arc length (to prevent singularity)

%% Possible Lane Change Intervals
% plot(lane_.l_inter(DSTART:DEND,1),'r-'); hold on; grid on;
% plot(lane_.r_inter(DSTART:DEND,1),'b-');
if sc == 1
    ub1 = find(lane_.state_idxs == DSTART + 1024);
    lb2 = find(lane_.state_idxs == DSTART + 1099);
    lane_.FactorValidIntvs = [1,ub1;lb2,length(imu_)];
    lane_.LC_dirs = {'right'};
elseif sc == 2
    lane_.FactorValidIntvs = [1,length(imu_)+1];
    lane.LC_dirs = {};
end
%% INS + GNSS Fusion 
% sol = struct();
% sol.basic = optimizer(imu_,gnss_,lane_,can_,snap,bias_,t_,covs_,'basic',options);
% sol.basic.optimize();
% 
% sol.basic.visualize();

%% INS + GNSS + WSS Fusion
sol = struct();
sol.partial = optimizer(imu_,gnss_,lane_,can_,snap,bias_,t_,covs_,'partial',options);
sol.partial.optimize();
%%
sol.partial.update('2-phase');


%% Processing VINS, DSO
% Read DSO
s3 = fopen('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vi_dso_result\scenario3.txt','r');
formatSpec = '%f %f %f %f %f %f %f %f';
Ssize = [8, Inf];
Rs3 = fscanf(s3,formatSpec,Ssize); Rs3 = Rs3';
DSO_S3 = [Rs3(:,2) Rs3(:,3)]; 

if sc == 1
    % Scenario 1: Curve
    vins_stereo_imu = readtable('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vins_fusion_stereo\scenario3_wIMU'); stfc = 1.4;
    vins_stereo = readtable('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vins_fusion_stereo\scenario3'); sfc = 1.75;
    vins_mono_imu = readtable('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vins_fusion_mono\scenario3'); mfc = 0.65;
    DSO = DSO_S3; dfc = 15; 
elseif sc == 2 
    % Scenario 2: Tunnel
    vins_stereo_imu = readtable('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vins_fusion_stereo\scenario4_wIMU'); stfc = 1.5;
    vins_stereo = readtable('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vins_fusion_stereo\scenario4'); sfc = 1;
    vins_mono_imu = readtable('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\vision_daejeon\vins_fusion_mono\scenario4'); mfc = 13.6;
    DSO = DSO_S3; dfc = 15; % fail
end



%% Post Processing VINS Data
vins = struct();
vins.stereoimu = struct(); vins.stereo = struct(); vins.monoimu = struct(); 

% VINS Stereo + IMU
vins.stereoimu.x = vins_stereo_imu.PosX;  vins.stereoimu.y = vins_stereo_imu.PosY; vins.stereoimu.heading = zeros(size(vins.stereoimu.x,1),1);
% scaling factor
d_vins_st_imu = sqrt((vins.stereoimu.x(2) - vins.stereoimu.x(1))^2 + (vins.stereoimu.y(2) - vins.stereoimu.y(1))^2);


for i=1:size(vins.stereoimu.x,1)-1
    vins.stereoimu.x(i+1) = lane.posx(DSTART) + stfc * (vins_stereo_imu.PosX(i+1) - vins.stereoimu.x(1)) * cos(lane.eul.z(DSTART)) - stfc * (vins_stereo_imu.PosY(i+1) - vins.stereoimu.y(1)) * sin(lane.eul.z(DSTART));
    vins.stereoimu.y(i+1) = lane.posy(DSTART) + stfc * (vins_stereo_imu.PosX(i+1) - vins.stereoimu.x(1)) * sin(lane.eul.z(DSTART)) + stfc * (vins_stereo_imu.PosY(i+1) - vins.stereoimu.y(1)) * cos(lane.eul.z(DSTART)); 
end

vins.stereoimu.x(1) = lane.posx(DSTART); vins.stereoimu.y(1) = lane.posy(DSTART);

% VINS Stereo
vins.stereo.x = vins_stereo.PosX; vins.stereo.y = vins_stereo.PosY; 

for i=1:size(vins.stereo.x,1)-1
    vins.stereo.x(i+1) = lane.posx(DSTART) + sfc * (vins_stereo.PosX(i+1) - vins.stereo.x(1)) * cos(lane.eul.z(DSTART)) - sfc * (vins_stereo.PosY(i+1) - vins.stereo.y(1)) * sin(lane.eul.z(DSTART));
    vins.stereo.y(i+1) = lane.posy(DSTART) + sfc* (vins_stereo.PosX(i+1) - vins.stereo.x(1)) * sin(lane.eul.z(DSTART)) + sfc * (vins_stereo.PosY(i+1) - vins.stereo.y(1)) * cos(lane.eul.z(DSTART));
end
vins.stereo.x(1) = lane.posx(DSTART); vins.stereo.y(1) = lane.posy(DSTART);
% VINS Mono + IMU
vins.monoimu.x = vins_mono_imu.PosX; vins.monoimu.y = vins_mono_imu.PosY;

for i=1:size(vins.monoimu.x,1)-1
    vins.monoimu.x(i+1) = lane.posx(DSTART) + mfc * (vins_mono_imu.PosX(i+1) - vins.monoimu.x(1)) * cos(lane.eul.z(DSTART)) - mfc * (vins_mono_imu.PosY(i+1) - vins.monoimu.y(1)) * sin(lane.eul.z(DSTART));
    vins.monoimu.y(i+1) = lane.posy(DSTART) + mfc * (vins_mono_imu.PosX(i+1) - vins.monoimu.x(1)) * sin(lane.eul.z(DSTART)) + mfc * (vins_mono_imu.PosY(i+1) - vins.monoimu.y(1)) * cos(lane.eul.z(DSTART));
end
vins.monoimu.x(1) = lane.posx(DSTART); vins.monoimu.y(1) = lane.posy(DSTART);

% DSO
dso = struct();
dso.x = DSO(:,1); dso.y = DSO(:,2);

for i=1:size(dso.x,1)-1
    dso.x(i+1) = lane.posx(DSTART) + dfc * (DSO(i+1,1) - dso.x(1)) * cos(lane.eul.z(DSTART) - pi/2) - dfc * (DSO(i+1,2) - dso.y(1)) * sin(lane.eul.z(DSTART) - pi/2);
    dso.y(i+1) = lane.posy(DSTART) + dfc * (DSO(i+1,1) - dso.x(1)) * sin(lane.eul.z(DSTART) - pi/2) + dfc * (DSO(i+1,2) - dso.y(1)) * cos(lane.eul.z(DSTART) - pi/2);
end
dso.x(1) = lane.posx(DSTART); dso.y(1) = lane.posy(DSTART);


%% Visualize with comparison 
if sc == 1
    sol.partial.visualizeRef(ref,vins.stereoimu,vins.stereo,vins.monoimu,dso);
elseif sc == 2
    sol.partial.visualizeRef(ref,vins.stereoimu,vins.stereo,vins.monoimu);
end



