classdef dataprocessor < handle
% DATAPROCESSOR - Read partially processed data by python and extract
% useful information or cluster data into usable format
% Processed data will be used in 'optimizer.m'
%
% Implemented by JinHwan Jeon, 2022
    properties
        raw_data = struct()
        data = struct()
        proc_data = struct() 
    end
    methods
        %% Constructor
        function obj = dataprocessor(varargin)
            imu = varargin{1};
            gnss = varargin{2};
            gnss.raw_idx = gnss.raw_idx + 1;
            can = varargin{3};
            lane = varargin{4};
            if length(varargin) == 5
                snap = varargin{5};
                obj.raw_data.snap = snap;
            end

            obj.raw_data.imu = imu;
            obj.raw_data.gnss = gnss;
            obj.raw_data.can = can;
            obj.raw_data.lane = lane;            
        end

        %% Process data
        function obj = process(obj)
            disp('[Beginning data processing...]')
            % Find non-static time intervals from wheel speed data
            valid_can_idxs = find(obj.raw_data.can.whl_spd > 0);
            
            obj.raw_data.can.valid_idxs = valid_can_idxs;
            obj.raw_data.can.valid_t = obj.raw_data.can.t(valid_can_idxs);
            
            lb = obj.raw_data.can.valid_idxs(1);
            obj.raw_data.t_intv = [];
            for i=1:length(obj.raw_data.can.valid_idxs)-1
                if obj.raw_data.can.valid_idxs(i+1) - obj.raw_data.can.valid_idxs(i) > 1                    
                    ub = obj.raw_data.can.valid_idxs(i);
                    obj.raw_data.t_intv = [obj.raw_data.t_intv; obj.raw_data.can.t(lb-1) obj.raw_data.can.t(ub+1)];
                    lb = obj.raw_data.can.valid_idxs(i+1);
                end
            end
            ub = obj.raw_data.can.valid_idxs(end);
            obj.raw_data.t_intv = [obj.raw_data.t_intv; obj.raw_data.can.t(lb-1) obj.raw_data.can.t(ub+1)];
            
            t_intv = obj.raw_data.t_intv;
            % Filter out invalid gnss measurements
            not_outage_idxs = find(obj.raw_data.gnss.hAcc < 5); % Remove GNSS measurements with hAcc larger than 3m
            obj.raw_data.gnss.t = obj.raw_data.gnss.t(not_outage_idxs);
            obj.raw_data.gnss.pos = obj.raw_data.gnss.pos(not_outage_idxs,:);
            obj.raw_data.gnss.hAcc = obj.raw_data.gnss.hAcc(not_outage_idxs);
            obj.raw_data.gnss.vAcc = obj.raw_data.gnss.vAcc(not_outage_idxs);
            obj.raw_data.gnss.vNED = obj.raw_data.gnss.vNED(not_outage_idxs,:);
            obj.raw_data.gnss.bearing = obj.raw_data.gnss.bearing(not_outage_idxs);
            obj.raw_data.gnss.raw_idx = obj.raw_data.gnss.raw_idx(not_outage_idxs);
            
            obj.data.gnss = {};
            obj.data.imu = {};
            
            gnss_ = struct();
            gnss_.t = [];
            gnss_.t_bias = 0;
            gnss_.lla_ref = obj.raw_data.gnss.lla_ref;
            gnss_.pos = [];
            gnss_.hAcc = [];
            gnss_.vAcc = [];
            gnss_.vNED = [];
            gnss_.bearing = [];
            gnss_.raw_idx = [];
            
            % Vehicle Local Frame: NED
            [accel_t,uniqueIdx] = unique(obj.raw_data.imu.accel.t');
            accel_meas = obj.raw_data.imu.accel.meas(uniqueIdx,:);
            accel_ = zeros(size(accel_meas));
            accel_(:,1) = -accel_meas(:,3);
            accel_(:,2) = accel_meas(:,2);
            accel_(:,3) = accel_meas(:,1); % need to consider gravitational effects
            
            [gyro_t,uniqueIdx] = unique(obj.raw_data.imu.gyro.t');
            gyro_meas = obj.raw_data.imu.gyro.meas(uniqueIdx,:);
            gyro_ = zeros(size(gyro_meas));
            gyro_(:,1) = -gyro_meas(:,3);
            gyro_(:,2) = gyro_meas(:,2);
            gyro_(:,3) = gyro_meas(:,1); 
            
%             disp('Displaying raw IMU data...')
%             figure(1);
%             plot(accel_(:,1),'r-'); hold on; grid on;
%             plot(accel_(:,2),'g-');
%             plot(accel_(:,3),'b-');
% 
%             figure(2);
%             plot(gyro_(:,1),'r-'); hold on; grid on;
%             plot(gyro_(:,2),'g-');
%             plot(gyro_(:,3),'b-');

            cnt = 1;
            
            % Need to fix if vehicle is not stopping near the final part of
            % dataset
            for i=1:length(obj.raw_data.gnss.t)
                if obj.raw_data.gnss.t(i) >= t_intv(cnt,1) && obj.raw_data.gnss.t(i) <= t_intv(cnt,2)
                    gnss_.t = [gnss_.t obj.raw_data.gnss.t(i)];
                    gnss_.pos = [gnss_.pos; obj.raw_data.gnss.pos(i,:)];
                    gnss_.hAcc = [gnss_.hAcc obj.raw_data.gnss.hAcc(i)];
                    gnss_.vAcc = [gnss_.vAcc obj.raw_data.gnss.vAcc(i)];
                    gnss_.vNED = [gnss_.vNED; obj.raw_data.gnss.vNED(i,:)];
                    gnss_.bearing = [gnss_.bearing obj.raw_data.gnss.bearing(i)];
                    gnss_.raw_idx = [gnss_.raw_idx obj.raw_data.gnss.raw_idx(i)];
                elseif obj.raw_data.gnss.t(i) > t_intv(cnt,2)
                    % Create IMU timestamps for interpolation
                    t_bias = gnss_.t_bias;
                    obj.data.gnss = [obj.data.gnss {gnss_}];

                    imu_ = struct();                    
                    t_lb = gnss_.t(1); t_ub = gnss_.t(end);
                    imu_.t = t_lb:0.01:t_ub;
                    imu_.t_bias = t_bias;
                    imu_.accel = interp1(accel_t,accel_,imu_.t,'linear','extrap');
                    imu_.gyro = interp1(gyro_t,gyro_,imu_.t,'linear','extrap');

                    obj.data.imu = [obj.data.imu {imu_}];
                    
                    if size(t_intv,1) == cnt
                        break;
                    else
                        
                        gnss_ = struct();
                        gnss_.t = [];
                        gnss_.t_bias = t_bias + t_intv(cnt+1,1) - t_intv(cnt,2);
                        gnss_.lla_ref = obj.raw_data.gnss.lla_ref;
                        gnss_.pos = [];
                        gnss_.hAcc = [];
                        gnss_.vAcc = [];
                        gnss_.vNED = [];
                        gnss_.bearing = [];
                        gnss_.raw_idx = [];
                        
                        cnt = cnt + 1;
                    end
                end
            end

            %% For each non-static interval, create additional timestamps for GNSS outages
            obj.data.base_unix_timestamp = obj.data.gnss{1}.t(1);
            base_t = obj.data.base_unix_timestamp;

            m = length(obj.data.gnss);
            obj.proc_data.full_t = [];
            obj.proc_data.t = [];
            obj.proc_data.gnss = struct();
%             obj.proc_data.gnss.t = obj.data.gnss{1}.t; % Change in the future for multi-static dataset
            obj.proc_data.gnss.t = [];
            obj.proc_data.gnss.state_idxs = [];
            obj.proc_data.gnss.pos = [];
            obj.proc_data.gnss.hAcc = [];
            obj.proc_data.gnss.vAcc = [];
            obj.proc_data.gnss.bearing = [];
            obj.proc_data.gnss.vNED = [];
            obj.proc_data.gnss.lla0 = obj.data.gnss{1}.lla_ref;

            cnt = 0;

            for i=1:m
                obj.proc_data.gnss.t = [obj.proc_data.gnss.t obj.data.gnss{i}.t - base_t - obj.data.gnss{i}.t_bias];
                obj.data.gnss{i}.full_t = [];
                n = length(obj.data.gnss{i}.raw_idx);
                t_bias = obj.data.gnss{i}.t_bias;
                
                % Augment GNSS data
                obj.proc_data.gnss.pos = [obj.proc_data.gnss.pos; obj.data.gnss{i}.pos];
                obj.proc_data.gnss.hAcc = [obj.proc_data.gnss.hAcc obj.data.gnss{i}.hAcc];
                obj.proc_data.gnss.vAcc = [obj.proc_data.gnss.vAcc obj.data.gnss{i}.vAcc];
                obj.proc_data.gnss.bearing = [obj.proc_data.gnss.bearing obj.data.gnss{i}.bearing];
                obj.proc_data.gnss.vNED = [obj.proc_data.gnss.vNED; obj.data.gnss{i}.vNED];

                for j = 1:n-1
                    % Add gnss idxs(state_relative)
                    obj.proc_data.gnss.state_idxs = [obj.proc_data.gnss.state_idxs cnt + length(obj.data.gnss{i}.full_t) + 1];

                    if obj.data.gnss{i}.raw_idx(j+1) - obj.data.gnss{i}.raw_idx(j) > 1 % If GNSS Outage is Detected
                        
                        % Error may occur if gnss timestamp spacing is
                        % exact, which is actually very unlikely.

                        obj.data.gnss{i}.full_t = [obj.data.gnss{i}.full_t obj.data.gnss{i}.t(j):0.1:obj.data.gnss{i}.t(j+1)-0.1];
                        obj.proc_data.t = [obj.proc_data.t (obj.data.gnss{i}.t(j):0.1:obj.data.gnss{i}.t(j+1)-0.1)-t_bias-base_t];
                        obj.proc_data.full_t = [obj.proc_data.full_t obj.data.gnss{i}.t(j):0.1:obj.data.gnss{i}.t(j+1)-0.1];
                    else
                        obj.data.gnss{i}.full_t = [obj.data.gnss{i}.full_t obj.data.gnss{i}.t(j)];
                        obj.proc_data.t = [obj.proc_data.t obj.data.gnss{i}.t(j) - t_bias - base_t];
                        obj.proc_data.full_t = [obj.proc_data.full_t obj.data.gnss{i}.t(j)];
                    end
                end
                obj.proc_data.gnss.state_idxs = [obj.proc_data.gnss.state_idxs cnt + length(obj.data.gnss{i}.full_t) + 1];
                obj.data.gnss{i}.full_t = [obj.data.gnss{i}.full_t obj.data.gnss{i}.t(end)];
                obj.proc_data.t = [obj.proc_data.t obj.data.gnss{i}.t(end) - t_bias - base_t];
                obj.proc_data.full_t = [obj.proc_data.full_t obj.data.gnss{i}.t(end)];

                cnt = cnt + length(obj.data.gnss{i}.full_t);
            end

            %% Cluster IMU measurements
            obj.proc_data.imu = {};

            for i=1:length(obj.data.gnss)
                full_t = obj.data.gnss{i}.full_t;
                imu_t = obj.data.imu{i}.t;
                t_bias = obj.data.gnss{i}.t_bias;
                
                accel_n = obj.data.imu{i}.accel;
                gyro_n = obj.data.imu{i}.gyro;
                
                if i ~= 1 % Continue adding from previous data fragment
                    if length(full_t) > 1
                        t_lb_ = full_t(1);
                        [~,idx_lb] = min(abs(imu_t - t_lb_));
                        t_ub = imu_t(idx_lb) - base_t - t_bias;
                        
                        imu_.t = t_lb:0.01:t_ub;
                        if t_ub - imu_.t(end) > 1e-6
                            imu_.t = [imu_.t t_ub];
                        end
                        
                        sampleT = [sampleT, t_ub];
                        sampleAccel = [sampleAccel; accel_n(idx_lb,:)];
                        sampleAccelX = sampleAccel(:,1);
                        sampleAccelY = sampleAccel(:,2);
                        sampleAccelZ = sampleAccel(:,3);
                        sampleGyro = [sampleGyro; gyro_n(idx_lb,:)];
                        sampleGyroX = sampleGyro(:,1);
                        sampleGyroY = sampleGyro(:,2);
                        sampleGyroZ = sampleGyro(:,3);
                        
                        interpedAccelX = interp1(sampleT,sampleAccelX,imu_.t);
                        interpedAccelY = interp1(sampleT,sampleAccelY,imu_.t);
                        interpedAccelZ = interp1(sampleT,sampleAccelZ,imu_.t);
                        interpedGyroX = interp1(sampleT,sampleGyroX,imu_.t);
                        interpedGyroY = interp1(sampleT,sampleGyroY,imu_.t);
                        interpedGyroZ = interp1(sampleT,sampleGyroZ,imu_.t);
                        
                        accel_N = [interpedAccelX(1:end-1)' interpedAccelY(1:end-1)' interpedAccelZ(1:end-1)'];
                        gyro_N = [interpedGyroX(1:end-1)' interpedGyroY(1:end-1)' interpedGyroZ(1:end-1)'];
                        
                        imu_.accel = accel_N;
                        imu_.gyro = gyro_N;
    
                        obj.proc_data.imu = [obj.proc_data.imu {imu_}];
    %                     disp(length(obj.proc_data.imu)) % Check if dt is well defined
                    end
                end

                for j=1:length(full_t)-1
                    t_lb = full_t(j); t_ub = full_t(j+1);
                    [~,idx_lb] = min(abs(imu_t - t_lb));
                    [~,idx_ub] = min(abs(imu_t - t_ub));
                    imu_ = struct();
                    imu_.t = [];
                    imu_.accel = [];
                    imu_.gyro = [];
                    
                    for k=idx_lb:idx_ub-1
                        imu_.t = [imu_.t imu_t(k) - base_t - t_bias];
                        imu_.accel = [imu_.accel; accel_n(k,:)];
                        imu_.gyro = [imu_.gyro; gyro_n(k,:)];
                    end
                    % Append last timestamp for computing dt
                    imu_.t = [imu_.t imu_t(idx_ub) - base_t - t_bias];
                    
                    obj.proc_data.imu = [obj.proc_data.imu {imu_}];
                end
                if length(full_t) == 1                    
                    n = length(obj.proc_data.imu);
                    obj.proc_data.full_t(n+1) = [];
                    obj.proc_data.t(n+1) = [];
                    rem_idxs = obj.proc_data.gnss.state_idxs(n+2:end);
                    obj.proc_data.gnss.state_idxs(n+1) = [];
                    obj.proc_data.gnss.state_idxs(n+1:end) = rem_idxs - 1;
                    obj.proc_data.gnss.t(n+1) = [];
                    obj.proc_data.gnss.pos(n+1,:) = [];
                    obj.proc_data.gnss.hAcc(n+1) = [];
                    obj.proc_data.gnss.vAcc(n+1) = [];
                    obj.proc_data.gnss.bearing(n+1) = [];
                    obj.proc_data.gnss.vNED(n+1,:) = [];
                else
                    % Between static time intervals
                    if i ~= length(obj.data.gnss)
                        
                        imu_ = struct();
                        imu_.t = [];
                        imu_.accel = [];
                        imu_.gyro = [];
                        
                        t_lb = imu_t(idx_ub) - base_t - t_bias;
                        sampleT = t_lb;
                        sampleAccel = accel_n(idx_ub,:);
                        sampleGyro = gyro_n(idx_ub,:);
    %                     for k=idx_ub:length(imu_t)
    %                         t_lb = 
    %                         imu_.t = [imu_.t imu_t(k) - base_t - t_bias];
    %                         imu_.accel = [imu_.accel; accel_n(k,:)];
    %                         imu_.gyro = [imu_.gyro; gyro_n(k,:)];
    %                     end
                    end
                end                
            end

            %% If there exists static interval, find bias values

            can_static_idxs = find(obj.raw_data.can.whl_spd == 0);
            can_static_timestamp = obj.raw_data.can.t(can_static_idxs);
            
            static_intvs = [];
            lb = can_static_timestamp(1);
            for i=1:length(can_static_idxs)-1
                if can_static_idxs(i+1) - can_static_idxs(i) > 1
                    ub = can_static_timestamp(i);
                    static_intvs = [static_intvs; lb ub];
                    lb = can_static_timestamp(i+1);
                end
            end
            static_intvs = [static_intvs; lb can_static_timestamp(end)];
            
            AccelBiasX = [];
            AccelBiasY = [];
            AccelBiasZ = [];
            GyroBiasX = [];
            GyroBiasY = [];
            GyroBiasZ = [];
            
            cnt = 1;
            for i=1:length(accel_t)
                if accel_t(i) > static_intvs(cnt,1) && accel_t(i) < static_intvs(cnt,2)
                    AccelBiasX = [AccelBiasX accel_(i,1)];
                    AccelBiasY = [AccelBiasY accel_(i,2)];
                    AccelBiasZ = [AccelBiasZ accel_(i,3) - 9.81]; % Change if data is in ENU format
                elseif accel_t(i) > static_intvs(cnt,2)
                    if cnt == size(static_intvs,1)
                        break;
                    else
                        cnt = cnt + 1;
                    end
                end
            end
            
            cnt = 1;
            for i=1:length(gyro_t)
                if gyro_t(i) > static_intvs(cnt,1) && gyro_t(i) < static_intvs(cnt,2)
                    GyroBiasX = [GyroBiasX gyro_(i,1)];
                    GyroBiasY = [GyroBiasY gyro_(i,2)];
                    GyroBiasZ = [GyroBiasZ gyro_(i,3)];
                elseif gyro_t(i) > static_intvs(cnt,2)
                    if cnt == size(static_intvs,1)
                        break;
                    else
                        cnt = cnt + 1;
                    end
                end
            end
            initAccelBias = [mean(AccelBiasX), mean(AccelBiasY), mean(AccelBiasZ)];
            initGyroBias = [mean(GyroBiasX), mean(GyroBiasY), mean(GyroBiasZ)];
            
            disp('Estimated Bias from static intervals')
            disp(['Accelerometer: ',num2str(initAccelBias)])
            disp(['Gyroscope: ',num2str(initGyroBias)])
            
            obj.proc_data.initBias = struct();
            obj.proc_data.initBias.accel = initAccelBias;
            obj.proc_data.initBias.gyro = initGyroBias;

            %% Find best match idx for CAN and Lane
            obj.proc_data.can = obj.raw_data.can;
            obj.proc_data.can.state_idxs = [];
            obj.proc_data.lane = obj.raw_data.lane;
            obj.proc_data.lane.state_idxs = [];
            for i=1:length(obj.proc_data.full_t)
                t = obj.proc_data.full_t(i);
                [~,idx_can] = min(abs(obj.raw_data.can.t - t));
                [~,idx_lane] = min(abs(obj.raw_data.lane.t - t));

                obj.proc_data.can.state_idxs = [obj.proc_data.can.state_idxs idx_can];
                obj.proc_data.lane.state_idxs = [obj.proc_data.lane.state_idxs idx_lane];
            end

            %% Process SNAP raw data
            % Perform mid-point integration to find timestamps for Nan
            % values

%             snap_t0 = obj.raw_data.snap.t(1);
%             [~,idx] = min(abs(obj.raw_data.can.t - snap_t0));
%             
%             vel = obj.raw_data.can.whl_spd;
%             t = obj.raw_data.can.t;
%             dist = 0;
%             
%             for i=idx:length(vel)-1
%                 dt = t(i+1) - t(i);
%                 avg_vel = 1/2 * (vel(i)+vel(i+1));
%                 dist = [dist dist(end) + avg_vel * dt];
%             end
%             obj.raw_data.snap.can_dist = dist;
% 
%             can_idxs = [];
%             est_t = [];
% 
%             for i=1:20
%                 L = obj.raw_data.snap.dist(i);
%                 [~,idx_] = min(abs(obj.raw_data.snap.can_dist - L));
%                 can_idxs = [can_idxs idx_];
%                 est_t = [est_t obj.raw_data.can.t(idx_ + idx - 1)];
%             end
% 
%             
%             obj.raw_data.snap.can_idxs = can_idxs;
%             
%             obj.raw_data.snap.est_t = est_t - snap_t0;
%             obj.raw_data.snap.ref_t = obj.raw_data.snap.t - snap_t0;

            
            %% Mark unstable lane data (lane change, etc)
            % Current rule: Do not perform optimization for lane data
            % during lane change
            % May not be using can information directly in the future
            % Lane Change Detection (윤진형)
            % 
            % * Need to differentiate lane change with left/right turns
            % How??
            %
            % Method 1: Lane Curvature
            % Method 2: IMU 


%             can_lane_idxs = [];
% 
%             for i=1:length(obj.raw_data.lane.t)
%                 lane_t = obj.raw_data.lane.t(i);
%                 [~,idx] = min(abs(lane_t - obj.raw_data.can.t));
%                 can_lane_idxs = [can_lane_idxs idx];
%             end
%             obj.raw_data.can.can_lane_idxs = can_lane_idxs;

            % Find Lane Change time intervals
            LC_intvs = [];
            LC_dirs = {};
            for i=2:length(obj.raw_data.can.t)
                % LC Start
                if (obj.raw_data.can.leftBlinker(i-1) == 0 && obj.raw_data.can.rightBlinker(i-1) == 0) && (obj.raw_data.can.leftBlinker(i) == 1 || obj.raw_data.can.rightBlinker(i) == 1)
                    t_lb = obj.raw_data.can.t(i);
                % LC End
                elseif (obj.raw_data.can.leftBlinker(i-1) == 1 || obj.raw_data.can.rightBlinker(i-1) == 1) && (obj.raw_data.can.leftBlinker(i) == 0 && obj.raw_data.can.rightBlinker(i) == 0)
                    t_ub = obj.raw_data.can.t(i);
                    LC_intvs = [LC_intvs; t_lb t_ub];

                    if obj.raw_data.can.leftBlinker(i-1) == 1
                        LC_dirs = [LC_dirs; {'left'}];
                    elseif obj.raw_data.can.rightBlinker(i-1) == 1
                        LC_dirs = [LC_dirs; {'right'}];
                    end
                end
            end

            % Scenario 1
%             LC_intvs = LC_intvs(1:end-1,:);
%             LC_dirs = LC_dirs(1:end-1);

            obj.proc_data.can.LC_intvs = LC_intvs;
            obj.proc_data.lane.LC_dirs = LC_dirs;
            

%             figure(1);
%             plot(obj.raw_data.can.t,obj.raw_data.can.leftBlinker,'r-'); hold on; grid on;
%             plot(obj.raw_data.can.t,obj.raw_data.can.rightBlinker,'b-');
%             for i=1:size(LC_intvs,1)
%                 T = LC_intvs(i,:);
%                 t = [T(1),T(1),T(2),T(2),T(1)];
%                 y = [0,1,1,0,0];
%                 plot(t,y,'k--',LineWidth=2);
%             end
            
            cnt = 1;
            laneFactorValidIdxs = true(1,length(obj.proc_data.full_t));
            for i=1:length(obj.proc_data.full_t)
                if obj.proc_data.full_t(i) >= LC_intvs(cnt,1) && obj.proc_data.full_t(i) <= LC_intvs(cnt,2)
                    laneFactorValidIdxs(i) = false;
                elseif obj.proc_data.full_t(i) > LC_intvs(cnt,2)
                    if cnt == size(LC_intvs,1)
                        break;
                    else
                        cnt = cnt + 1;
                    end                    
                end
            end
            
            % laneFactorValidIntvs: Valid state idxs for lane factor 
            % Other than these indices, lane variables are not defined
            % Format: [lower_bound_index, upper_bound_index, number of valid indices in this interval]
            % Use 'sum.m' to find total number of valid indices when
            % computing the height of Block Jacobian or Measurement
            % Residual.

            lb = 1;
            laneFactorValidIntvs = [];
            for i=2:length(laneFactorValidIdxs)
                if laneFactorValidIdxs(i-1) == 1 && laneFactorValidIdxs(i) == 0
                    laneFactorValidIntvs = [laneFactorValidIntvs; lb i-1 i-1-lb+1];
                elseif laneFactorValidIdxs(i-1) == 0 && laneFactorValidIdxs(i) == 1
                    lb = i;
                end
            end

            laneFactorValidIntvs = [laneFactorValidIntvs; lb length(laneFactorValidIdxs) length(laneFactorValidIdxs)-lb+1];
            
            obj.proc_data.lane.FactorValidIntvs = laneFactorValidIntvs;
            
            % Augmented Indices
            org = 1:length(laneFactorValidIdxs);
            obj.proc_data.lane.FactorValidIdxs = org(laneFactorValidIdxs);

            %% Remove Pure Left/Right Turns from lane change
             % Remove last "right turn" as it is not a lane change scenario
            % This part should be changed for other dataset



            
            % Scenario 2
            laneFactorValidIntvs(4,:) = [laneFactorValidIntvs(4,1), laneFactorValidIntvs(5,2), ...
                                         laneFactorValidIntvs(5,2) - laneFactorValidIntvs(4,1)+1];
            laneFactorValidIntvs(5,:) = [];
            laneFactorValidIntvs(end,1) = laneFactorValidIntvs(end,1) + 15;
            obj.proc_data.lane.FactorValidIntvs = laneFactorValidIntvs;
            obj.proc_data.lane.LC_dirs(4) = [];

        end
        
        %% Visualize Processed data
        function obj = visualize(obj)
            t = [];
            accel = [];
            gyro = [];

            for i=1:length(obj.data.imu)
                t = [t obj.data.imu{i}.t];
                accel = [accel; obj.data.imu{i}.accel];
                gyro = [gyro; obj.data.imu{i}.gyro];
            end
    
            figure(1); hold on; grid on; axis tight;
            plot(t',accel(:,1),'r-');
            plot(t',accel(:,2),'g-');
            plot(t',accel(:,3),'b-');
            xlabel('Timestamp'); ylabel('Acceleration (m/s^2)');
            title('Interpolated Acceleration Data')
            legend('Ax','Ay','Az')

            figure(2); hold on; grid on; axis tight;
            plot(t',gyro(:,1),'r-');
            plot(t',gyro(:,2),'g-');
            plot(t',gyro(:,3),'b-');
            xlabel('Timestamp'); ylabel('Gyro (rad/s)');
            title('Interpolated Gyroscope Data')
            legend('Wx','Wy','Wz')
        end

    end
end