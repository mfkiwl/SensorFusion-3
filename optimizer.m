classdef optimizer < handle
% OPTIMIZER - Sensor Fusion Optimizer Module
%   All input data are from the dataprocessor module
%   [Usage Format: an example]
%   sol = OPTIMIZER(imu,gnss,lane,can,imu_bias,t,prior_covs,'basic',options)
%
%   [imu, gnss, lane, can, imu_bias, t] : Pre-processed data
% 
%   [Mode] : 'basic', 'partial', 'full', '2-phase'
%   * Basic: INS + GNSS Fusion
%   * Partial: INS + GNSS + WSS Fusion
%   * Full, 2-phase: INS + GNSS + WSS + Lane Detection Fusion 
%   - Currently, 'Full' mode is not implemented due to stability issues
% 
%   [Options] : struct variable containing NLS optimization options
%   * options.CostThres: Cost difference threshold
%   * options.StepThres: Step size threshold
%   * options.IterThres: Iteration limit
%   * options.Algorithm: GN(Gauss-Newton),LM(Levenberg-Marquardt),TR(Trust-Region)
%   ** GN and LM : no guarantee of convergence to local minima
%   Recommended to use GN or TR 
%   When using TR as solver algorithm, need to define parameters.
%   -TR parameters at 'main.m' should work fine
% 
%   [Methods] : using OPTIMIZER 
%   * optimize() : Find optimized solution to SNLS problem
%   * visualize() : Plot optimization results
%   * update() : Update optimization mode and use previous optimized results
%   * Other function methods are not designed to be used outside of this script 
%
%   Usage example is introduced in "main.m". 
%
%   Implemented by JinHwan Jeon, 2022

    properties (Access = public)
        imu
        gnss
        lane
        can
        snap
        imu_bias
        t
        covs
        mode
        params % vehicle parameters (for partial or above modes)
        % * L: Relative vector from vehicle rear axle to IMU in (body frame)

        states = {} % state variables saver
        % * R: Rotation Matrix from vehicle body frame to world frame
        % * V: Velocity vector in world local frame
        % * P: Position vector in world local frame
        % ----- For 'partial' or above modes -----
        % * wsf: Wheel Scaling Factor(Ratio for considering effective radius)
    
        bias = {} % bias variables saver
        bgd_thres = 1e-2; % Gyro bias repropagation threshold
        bad_thres = 1e-1; % Accel bias repropagation threshold

        map % map variable for 'full' or '2-phase' mode
        phase % phase number for 
        opt = struct() % Optimized results
    end
    
    %% Public Methods
    methods (Access = public) 
        %% Constructor
        function obj = optimizer(varargin)
            obj.imu = varargin{1};
            obj.gnss = varargin{2};
            obj.lane = varargin{3};
            obj.can = varargin{4};
            obj.snap = varargin{5};
            obj.imu_bias = varargin{6};
            obj.t = varargin{7};
            obj.covs = varargin{8};
            obj.mode = varargin{9};            
            obj.opt.options = varargin{10};

            % Read preview number for full, 2-phase optimization
            
            n = length(obj.imu)+1;
            for i=1:n
                bias = struct();
                bias.bg = obj.imu_bias.gyro';
                bias.ba = obj.imu_bias.accel';
                bias.bgd = zeros(3,1);
                bias.bad = zeros(3,1);

                obj.bias = [obj.bias {bias}];
            end
            
            % Create Initial value for vehicle parameters
            if ~strcmp(obj.mode,'basic')
                obj.params =  [1.5; 0; 1];  
                obj.phase = 1;
            end

            % Perform Pre-integration
            obj.integrate(1:n-1);

            % Perform INS propagation
            obj.ins();

            % If mode is full, 2-phase, create initial value for lane
            % variables
        end
        
        %% Update Optimization Mode (for 2-phase)
        function obj = update(obj,str)
            if ismember(str,{'basic','partial','full','2-phase'})
                obj.mode = str;
            else
                error('Mode Selection: basic, partial, full, 2-phase')
            end
            
            if strcmp(obj.mode,'2-phase')
                disp("[Updating mode to 2-phase, creating Map information...]")
                % Create and initialize Map
                % Initial segmentation and data association is done in this
                % step automatically
                % Complete arc fitting is done in this step
                obj.map = ArcMap(obj.states,obj.lane,obj.lane.prob_thres,obj.opt.LeftCov,obj.opt.RightCov);
            end           
        end

        %% Optimize
        function obj = optimize(obj)
            
            % Solve Sparse Nonlinear Least Squares in Batch Manner                     

            disp('[Optimization Starts...]')

            n = length(obj.states);

            if ~strcmp(obj.mode,'2-phase') && ~strcmp(obj.mode,'full')
                % Run optimization depending on algorithm options
                %
                % Phase 1 Optimization
                %
                % When lane data is not used, only phase 1 optimization is
                % performed. 
                % [Basic Mode]: INS + GNSS Fusion
                % [Partial Mode]: INS + GNSS + WSS Fusion

                if strcmp(obj.mode,'basic')
                    obj.opt.x0 = zeros(15*n,1);
                elseif strcmp(obj.mode,'partial')
                    obj.opt.x0 = zeros(16*n,1);                            
                end

                if strcmp(obj.opt.options.Algorithm,'GN')
                    obj.GaussNewton();
                elseif strcmp(obj.opt.options.Algorithm,'LM')
                    obj.LevenbergMarquardt();
                elseif strcmp(obj.opt.options.Algorithm,'TR')
                    obj.TrustRegion();
                end
            else
                % Lane data is used for optimization
                % Current model: Continuous Piecewise Arc Spline Lane Model
                % 
                % Phase 2 Optimization
                % 
                % =================<Optimization step>=================
                % 1. Perform data association with current optimized
                % vehicle states and arc parameters
                %
                % 2. Fully converge the problem with iterative data
                % association
                %
                % 3. Validate optimization results and add subsegments if
                % needed. 
                %
                % 4. Until all segments are valid, repeat steps 1 ~ 3.
                % =====================================================
                %
                
                obj.phase = 2;
                obj.opt.options.Algorithm = 'TR';
                obj.opt.options.IterThres = inf; % Change if needed
                obj.opt.options.CostThres = 1;
                valid = false;

                while ~valid
                    
                    obj.updateTracker();                    
                    obj.opt.x0 = zeros(16*n + obj.opt.seg_tracker(end),1);

                    % Step 2: Full Convergence 
                    if strcmp(obj.opt.options.Algorithm,'GN')
                        obj.GaussNewton();
                    elseif strcmp(obj.opt.options.Algorithm,'LM')
                        obj.LevenbergMarquardt();
                    elseif strcmp(obj.opt.options.Algorithm,'TR')
                        obj.TrustRegion();
                    end                    
                    
                    % Step 3: Validate and perform data association 
                    % Create new arc segments if needed 
                    
                    obj.map.validate();
                    
                    % Step 1 is done inside the validation function
                    valid = obj.map.valid_flag;

                    % Step 4: Repeat 1 ~ 3 until all segments are valid
                    % Show optimization results after optimization
                    % 
                    % Number of Sub-segments for each segments
                    % Summary for each Sub-segments: Curvature, Arc Length
                    % Summary for each Segments: Curvature, Arc Length
                end               
            end
        end
         
        %% Visualize
        function visualize(obj)
            n = length(obj.states);
            
            R = []; V = []; P = []; Bg = []; Ba = [];
            V_b = []; whl_spd = []; S = [];
            left = {}; right = {};

            for i=1:n
                R_ = obj.states{i}.R;
                V_ = obj.states{i}.V;
                P_ = obj.states{i}.P;
                
                Left_ = []; Right_ = [];
                % Transform lane points to world frame
                
                 for j=1:obj.lane.prev_num
                    left_ = obj.states{i}.left(:,j);
                    right_ = obj.states{i}.right(:,j);

                    Left_ = [Left_ P_ + R_ * left_];
                    Right_ = [Right_ P_ + R_ * right_];

                    
                end
                left = [left {Left_}];
                right = [right {Right_}];

%                 if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
%                     for j=1:obj.lane.prev_num
%                         left_ = obj.states{i}.left(:,j);
%                         right_ = obj.states{i}.right(:,j);
%     
%                         Left_ = [Left_ P_ + R_ * left_];
%                         Right_ = [Right_ P_ + R_ * right_];
%     
%                         
%                     end
%                     left = [left {Left_}];
%                     right = [right {Right_}];
%                 end
                
                R = [R dcm2rpy(R_)];
                V = [V V_]; P = [P P_];
                V_b = [V_b R_' * V_];
                Bg_ = obj.bias{i}.bg;
                Ba_ = obj.bias{i}.ba;
                Bg = [Bg Bg_]; Ba = [Ba Ba_];

                can_idx = obj.can.state_idxs(i);
                whl_spd = [whl_spd obj.can.whl_spd(can_idx)];
                if ~strcmp(obj.mode,'basic')
                    S = [S obj.states{i}.WSF];
                end
            end
            gnss_pos = lla2enu(obj.gnss.pos,obj.gnss.lla0,'ellipsoid');
            

            figure(1); hold on; grid on; axis equal;
%             p_est = plot3(P(1,:),P(2,:),P(3,:),'r.');
%             p_gnss = plot3(gnss_pos(:,1),gnss_pos(:,2),gnss_pos(:,3),'b.');
%             
%             for i=1:n
%                 left_ = left{i}; right_ = right{i};
%                 plot3(left_(1,:),left_(2,:),left_(3,:));
%                 plot3(right_(1,:),right_(2,:),right_(3,:));
%             end
            p_est = plot(P(1,:),P(2,:),'r.');
            p_gnss = plot(gnss_pos(:,1),gnss_pos(:,2),'b.');
            
            for i=1:n
                left_ = left{i}; right_ = right{i};
                plot(left_(1,1),left_(2,1),'cx');
                plot(right_(1,1),right_(2,1),'mx');
            end

            for i=1:size(obj.lane.FactorValidIntvs,1)-1
                ub = obj.lane.FactorValidIntvs(i,2);
                plot(obj.states{ub}.P(1),obj.states{ub}.P(2),'g+');
            end

%             if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
%                 for i=1:n
%                     left_ = left{i}; right_ = right{i};
%                     plot(left_(1,:),left_(2,:));
%                     plot(right_(1,:),right_(2,:));
%                 end
%             end

            xlabel('Global X'); ylabel('Global Y');
            title('Optimized Trajectory 3D Comparison')
            legend([p_est,p_gnss],'Estimated Trajectory','GNSS Measurements')
            
            figure(2); 
%             snap_lla = [obj.snap.lat, obj.snap.lon, zeros(size(obj.snap.lat,1),1)];                       
%             p_snap = geoplot(dataset.raw_data.snap.lat,dataset.raw_data.snap.lon,'r.');
            p_gnss = geoplot(obj.gnss.pos(:,1),obj.gnss.pos(:,2),'b.'); grid on; hold on;
            P_lla = enu2lla(P',obj.gnss.lla0,'ellipsoid');                        
            obj.opt.P_lla = P_lla;
            p_est = geoplot(P_lla(:,1),P_lla(:,2),'r.'); 
            
            if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
                n = length(obj.map.arc_segments);
                for i=1:n
                    m = length(obj.map.arc_segments{i}.kappa);
                    seg = obj.map.arc_segments{i};
                    heading = seg.tau0;
                    SegPoints = [seg.x0;
                                 seg.y0];
                    for j=1:m
                        kappa = seg.kappa(j); L = seg.L(j);
                        headingPrev = heading;
                        heading = heading + kappa * L;
                        headingCurr = heading;
    
                        heading_ = linspace(headingPrev,headingCurr,1e3);
                        addedSegPoints = SegPoints(:,end) + 1/kappa * [sin(heading_) - sin(headingPrev);
                                                                       -cos(heading_) + cos(headingPrev)];
                        SegPoints = [SegPoints addedSegPoints];
                        
                        lla_node1 = enu2lla([addedSegPoints(1,1),addedSegPoints(2,1),0],obj.gnss.lla0,'ellipsoid');
                        lla_node2 = enu2lla([addedSegPoints(1,end),addedSegPoints(2,end),0],obj.gnss.lla0,'ellipsoid');
                        p_node = geoplot(lla_node1(1),lla_node1(2),'co');
                        geoplot(lla_node2(1),lla_node2(2),'co');
                    end
                    
                    lla_seg = enu2lla([SegPoints' zeros(size(SegPoints,2),1)],obj.gnss.lla0,'ellipsoid');
                    lla_node3 = enu2lla([SegPoints(1,1),SegPoints(2,1),0],obj.gnss.lla0,'ellipsoid');
                    lla_node4 = enu2lla([SegPoints(1,end),SegPoints(2,end),0],obj.gnss.lla0,'ellipsoid');
    
                    p_lane = geoplot(lla_seg(:,1),lla_seg(:,2),'k-');
                    p_seg = geoplot(lla_node3(1),lla_node3(2),'ms');
                    geoplot(lla_node4(1),lla_node4(2),'ms');
                    
                    
                end
                geobasemap satellite

                title('Optimized Trajectory Comparison')
                legend([p_est,p_gnss,p_lane,p_seg,p_node],...
                       'Optimized Vehicle Trajectory','GNSS Measurements',...
                       'Continuous Piecewise Arc Spline',...
                       'Arc Segment Node','Arc Sub-Segment Node')
%                 title('Scenario Trajectory')
%                 legend([p_gnss],'GNSS Measurements')
            else
                geobasemap satellite

                title('Optimized Trajectory Comparison')
                legend([p_est,p_gnss],...
                       'Optimized Vehicle Trajectory','GNSS Measurements')
            end

%             geoplot(obj.snap.lat,obj.snap.lon,'c.','MarkerSize',8)

            

            figure(3); hold on; grid on;
            plot(obj.t,V(1,:),'r-');
            plot(obj.t,V(2,:),'g-');
            plot(obj.t,V(3,:),'b-');
            xlabel('Time(s)'); ylabel('V(m/s)')
            title('Optimized Velocity')
            legend('Vx','Vy','Vz')
            
            figure(4); hold on; grid on;
            plot(obj.t,V_b(1,:),'r-');
            plot(obj.t,V_b(2,:),'g-');
            plot(obj.t,V_b(3,:),'b-');
            plot(obj.t,whl_spd,'c-');
            xlabel('Time(s)'); ylabel('V(m/s)')
            title('Optimized Velocity (Body Frame)')
            legend('Vx_{b}','Vy_{b}','Vz_{b}','Wheel Speed','Location','NorthWest')

            figure(5); hold on; grid on;

            psi = wrapToPi(deg2rad(90 - obj.gnss.bearing));

            plot(obj.t,R(1,:),'r-');
            plot(obj.t,R(2,:),'g-');
            plot(obj.t,R(3,:),'b-');
            plot(obj.gnss.t-obj.gnss.t(1),psi,'c.')
            xlabel('Time(s)'); ylabel('Angle(rad)')
            title('Optimized RPY Angle Comparison')
            legend('Roll','Pitch','Yaw','GNSS Yaw')
            
            figure(6); hold on; grid on;
            plot(obj.t,Bg(1,:),'r-');
            plot(obj.t,Bg(2,:),'g-');
            plot(obj.t,Bg(3,:),'b-');
            xlabel('Time(s)'); ylabel('Gyro Bias(rad/s)')
            title('Optimized Gyro Bias Comparison')
            legend('wbx','wby','wbz')

            figure(7); hold on; grid on;
            plot(obj.t,Ba(1,:),'r-');
            plot(obj.t,Ba(2,:),'g-');
            plot(obj.t,Ba(3,:),'b-');
            xlabel('Time(s)'); ylabel('Accel Bias(m/s^2)')
            title('Optimized Accel Bias Comparison')
            legend('abx','aby','abz')

            if ~strcmp(obj.mode,'basic')
                figure(8); 
                plot(obj.t,S); grid on;
                xlabel('Time(s)'); ylabel('Wheel Speed Scaling Factor')
                title('Optimized Wheel Speed Scaling Factor')
            end
        end
        
        %% Visualize with HD Map
        function visualizeHD(obj,Lane)
            figure(20); hold on; grid on; axis equal;
            p_l = mapshow(Lane);             

            % Convert vehicle states to UTM-K coords            
            
            ENU = [];
            for i=1:length(obj.states)
                ENU = [ENU; obj.states{i}.P'];                
            end
            LLA = enu2lla(ENU,obj.gnss.lla0,'ellipsoid');
            p = projcrs(32652); % Projected Coordinate Reference System
            [x,y] = projfwd(p,LLA(:,1),LLA(:,2));
            p_v = plot(x,y,'k.');

            xlabel('UTM-K X (m)');ylabel('UTM-K Y (m)');
            title('HD Map')
            legend([p_l, p_v],'HD Map Lanes','Vehicle Trajectory')
        end        
        
        %% Plot lane point std values
        function plotConfEllipse(obj,p)
            % 0 < p < 1 : Confidence Level
            % Left Lane
            figure(1); grid on; axis equal; hold on;
            % i=1:length(obj.states)
            for i=700:1000
                R = obj.states{i}.R;
                P = obj.states{i}.P;
                Left = obj.states{i}.left;
                Right = obj.states{i}.right;
                % 1:obj.lane.prev_num
                for j=1:3
                    LeftL = P + R * Left(:,j);
                    RightL = P + R * Right(:,j);

                    CovL = reshape(obj.opt.LeftCov(4*i-3:4*i,j),2,2);
                    CovR = reshape(obj.opt.RightCov(4*i-3:4*i,j),2,2);

                    PC_L = obj.ConfidenceEllipse(LeftL(1:2),CovL,p);
                    PC_R = obj.ConfidenceEllipse(RightL(1:2),CovR,p);

                    p_l = plot(LeftL(1),LeftL(2),'kx');
                    plot(RightL(1),RightL(2),'kx');
                    p_cl = plot(PC_L(1,:),PC_L(2,:),'m--');
                    p_cr = plot(PC_R(1,:),PC_R(2,:),'c--');
                end
            end
            
            xlabel('Global X'); ylabel('Global Y');
            title(['Lane Point ',num2str(100 * p,2),'% Confidence Ellipse'])
            legend([p_l,p_cl,p_cr], ... 
                   'Lane Points','Left Lane Confidence Ellipse','Right Lane Confidence Ellipse')
            
        end

        %% Plot Vehicle position confidence ellipse
        function plotConfEllipseV(obj,p)
            figure(1); grid on; axis equal; hold on;
            for i=1:length(obj.states)
                P = obj.states{i}.P;
                xycov = full(obj.opt.cov(9*(i-1)+7:9*(i-1)+8,9*(i-1)+7:9*(i-1)+8));

                pc = obj.ConfidenceEllipse(P(1:2),xycov,p);
                plot(P(1),P(2),'k.');
                plot(pc(1,:),pc(2,:),'c--');
            end
        end

    end
    
    %% Private Methods
    methods (Access = private)
        %% Pre-integration using IMU Measurements, given IMU idxs
        function obj = integrate(obj,idxs)
%             disp(length(idxs))
            m = length(idxs);
            nbg_cov = obj.covs.imu.GyroscopeNoise;
            nba_cov = obj.covs.imu.AccelerometerNoise;
            n_cov = blkdiag(nbg_cov,nba_cov);

            for i=1:m
                % Integrate for each IMU clusters
                idx = idxs(i);
                n = size(obj.imu{idx}.accel,1);
                t_ = obj.imu{idx}.t;
                a = obj.imu{idx}.accel';
                w = obj.imu{idx}.gyro';
                ab = obj.bias{idx}.ba;
                wb = obj.bias{idx}.bg;

                delRik = eye(3); delVik = zeros(3,1); delPik = zeros(3,1);
                JdelRik_bg = zeros(3,3);
                JdelVik_bg = zeros(3,3);
                JdelVik_ba = zeros(3,3);
                JdelPik_bg = zeros(3,3);
                JdelPik_ba = zeros(3,3);

                noise_cov = zeros(9,9);
                noise_vec = zeros(9,1);

                dtij = 0;

                for k=1:n
                    dt_k = t_(k+1) - t_(k);
                    dtij = dtij + dt_k;

%                     n_k = mvnrnd(zeros(6,1),1/dt_k * n_cov)';
                    a_k = a(:,k) - ab;
                    w_k = w(:,k) - wb;

                    delRkkp1 = Exp_map(w_k * dt_k);
                    [Ak, Bk] = obj.getCoeff(delRik,delRkkp1,a_k,w_k,dt_k);

%                     noise_vec = Ak * noise_vec + Bk * n_k;
                    noise_cov = Ak * noise_cov * Ak' + Bk * (1/dt_k * n_cov) * Bk';
                    
                    % IMU measurement propagation
                    delPik = delPik + delVik * dt_k + 1/2 * delRik * a_k * dt_k^2;
                    delVik = delVik + delRik * a_k * dt_k;
                    
                    % Jacobian Propagation
                    JdelPik_ba = JdelPik_ba + JdelVik_ba * dt_k - 1/2 * delRik * dt_k^2;
                    JdelPik_bg = JdelPik_bg + JdelVik_bg * dt_k - 1/2 * delRik * skew(a_k) * JdelRik_bg * dt_k^2;
                    JdelVik_ba = JdelVik_ba - delRik * dt_k;
                    JdelVik_bg = JdelVik_bg - delRik * skew(a_k) * JdelRik_bg * dt_k;
                    JdelRik_bg = delRkkp1' * JdelRik_bg - RightJac(w_k * dt_k) * dt_k;
                    
                    delRik = delRik * Exp_map(w_k * dt_k);
                    
                end

                obj.imu{idx}.JdelRij_bg = JdelRik_bg;
                obj.imu{idx}.JdelVij_bg = JdelVik_bg;
                obj.imu{idx}.JdelVij_ba = JdelVik_ba;
                obj.imu{idx}.JdelPij_bg = JdelPik_bg;
                obj.imu{idx}.JdelPij_ba = JdelPik_ba;

                obj.imu{idx}.delRij = delRik;
                obj.imu{idx}.delVij = delVik;
                obj.imu{idx}.delPij = delPik;
                obj.imu{idx}.Covij = noise_cov;
                obj.imu{idx}.nij = noise_vec;

                lambdas = eig(noise_cov);
                if ~isempty(find(lambdas <= 0,1))
                    disp(idx)
                    disp(lambdas)
                    error('Computed IMU Covariance is not positive definite')
                end
                
%                 obj.imu{idx}.delRijn = delRik * Exp_map(-noise_vec(1:3));
%                 obj.imu{idx}.delVijn = delVik - noise_vec(4:6);
%                 obj.imu{idx}.delPijn = delPik - noise_vec(7:9);

                obj.imu{idx}.dtij = dtij;

            end
        end
            
        %% INS Propagation 
        function obj = ins(obj)
        % INS Propagation in the world frame using pre-integrated values

        % * Vehicle Body Frame: [x, y, z] = [Forward, Left, Up]
        % * World Frame: [x, y, z] = [East, North, Up]
        % * Camera Frame(Phone): https://source.android.com/docs/core/sensors/sensor-types
        % Note that IMU (Android) measurements have gravitational effects,
        % which should be considered when modeling vehicle 3D kinematics 

        % # Vehicle State Variables
        % # R: Body-to-world frame rotational matrix
        % # V: World frame velocity vector
        % # P: World frame position vector (local frame coords, not geodetic)
            disp('[INS Propagation...]')
            th = pi/2 - pi/180 * obj.gnss.bearing(1);
            R = [cos(th) -sin(th) 0;
                 sin(th) cos(th)  0;
                 0       0        1];
            Vned = obj.gnss.vNED(1,:);
            V = [Vned(2); Vned(1); -Vned(3)];
            P = lla2enu(obj.gnss.pos(1,:),obj.gnss.lla0,'ellipsoid')';
            
            grav = [0;0;-9.81];

            state_ = struct();
            state_.R = R; state_.V = V; state_.P = P;
            if ~strcmp(obj.mode,'basic')
                state_.WSF = 1;
            end
            lane_idx = obj.lane.state_idxs(1);
            prev_num = obj.lane.prev_num;

            state_.left = [0:10:10*(prev_num-1);
                           obj.lane.ly(lane_idx,1:prev_num);
                           obj.lane.lz(lane_idx,1:prev_num)];

            state_.right = [0:10:10*(prev_num-1);
                            obj.lane.ry(lane_idx,1:prev_num);
                            obj.lane.rz(lane_idx,1:prev_num)];

            obj.states = [obj.states {state_}];
            
            % For initial state, add bias information
            state_.bg = obj.bias{1}.bg;
            state_.bgd = obj.bias{1}.bgd;
            state_.ba = obj.bias{1}.ba;
            state_.bad = obj.bias{1}.bad;
            if ~strcmp(obj.mode,'basic')
                state_.L = obj.params;
            end

            obj.opt.init_state = state_;
            

            n = length(obj.imu);

            for i=1:n
                delRij = obj.imu{i}.delRij;
                delVij = obj.imu{i}.delVij;
                delPij = obj.imu{i}.delPij;
                dtij = obj.imu{i}.dtij;
                
                P = P + V * dtij + 1/2 * grav * dtij^2 + R * delPij;
                V = V + grav * dtij + R * delVij;
                R = R * delRij;

                state_ = struct();
                state_.R = R; state_.V = V; state_.P = P;
                
                if ~strcmp(obj.mode,'basic')
                    state_.WSF = 1;                    
                end
                lane_idx = obj.lane.state_idxs(i+1);  

                % Perhaps change this part so that lane states are only
                % created for only valid intvs
                % (Exclude Lane Changing intervals)
                state_.left = [0:10:10*(prev_num-1);
                               obj.lane.ly(lane_idx,1:prev_num);
                               obj.lane.lz(lane_idx,1:prev_num)];

                state_.right = [0:10:10*(prev_num-1);
                                obj.lane.ry(lane_idx,1:prev_num);
                                obj.lane.rz(lane_idx,1:prev_num)];

                obj.states = [obj.states {state_}];
            end

        end

        %% Gauss-Newton Method
        function obj = GaussNewton(obj)
            % Gauss-Newton method shows very fast convergence compared to 
            % Trust-Region method, but suffers from low convergence stability
            % near the local minima (oscillation frequently observed)
            % Oscillation detection was implemented to reduce meaningless
            % iterations. 

            disp('[SNLS solver: Gauss-Newton Method]')
            disp(['Current Phase Number: ',num2str(obj.phase)])
            fprintf(' Iteration     f(x)           step\n');
            formatstr = ' %5.0f   %13.6g  %13.6g ';
            x0 = obj.opt.x0;
            [res,jac] = obj.cost_func(x0);
            prev_cost = res' * res;
            str = sprintf(formatstr,0,prev_cost,norm(x0));
            disp(str)

            i = 1;
            
            cost_stack = prev_cost; % Stack costs for oscillation detection

            while true
                A = jac' * jac;
                b = -jac' * res;

                x0 = A \ b; 
                                
                [res,jac] = obj.cost_func(x0);

                cost = res' * res;
                cost_stack = [cost_stack cost];
                step_size = norm(x0);
                
                str = sprintf(formatstr,i,cost,step_size);
                disp(str)
                
                % Ending Criterion
                obj.opt.flags = [];
                obj.opt.flags = [obj.opt.flags abs(prev_cost - cost) > obj.opt.options.CostThres];
                obj.opt.flags = [obj.opt.flags step_size > obj.opt.options.StepThres];
                obj.opt.flags = [obj.opt.flags i < obj.opt.options.IterThres];
                
                % Check for oscillation around the local minima
                if length(cost_stack) >= 5
                    osc_flag = obj.DetectOsc(cost_stack);
                else
                    osc_flag = false;
                end

                if length(find(obj.opt.flags)) ~= length(obj.opt.flags) % If any of the criterion is not met, end loop
                    
                    obj.retract(x0,'final');

                    disp('[Optimization Finished...]')
                    idx = find(~obj.opt.flags,1);
                    
                    if idx == 1
                        disp(['Current cost difference ',num2str(abs(prev_cost-cost)),' is below threshold: ',num2str(obj.opt.options.CostThres)])
                    elseif idx == 2
                        disp(['Current step size ',num2str(step_size),' is below threshold: ',num2str(obj.opt.options.StepThres)])
                    elseif idx == 3
                        disp(['Current iteration number ',num2str(i),' is above threshold: ',num2str(obj.opt.options.IterThres)])
                    end

                    break;
                elseif osc_flag
                    obj.retract(x0,'final');

                    disp('[Optimization Finished...]')
                    disp('Oscillation about the local minima detected')
                    break;
                else
                    i = i + 1;
                    prev_cost = cost;
                end
                
            end

            % Sparse Inverse to extract lane covariance
            obj.opt.info = jac' * jac;
            disp('[Computing Sparse Inverse of Information Matrix...]')
            obj.opt.cov = sparseinv(obj.opt.info);
            obj.extractLaneCov();
        end

        %% Levenberg-Marquardt Method
        function obj = LevenbergMarquardt(obj)
            % Levenberg-Marquardt algorithm
            % Implemented "The Levenberg-Marquardt algorithm for nonlinear least squares curve-fitting problems"
            % by Henri P. Gavin
            % 
            % Implemented in MATLAB by JinHwan Jeon, 2022
            %
            % Original paper provides 3 possible update methods for L-M 
            % parameter but only the first update method is implemented
            % Levenberg-Marquardt Algorithm seems to be not an appropriate
            % solver for large sparse optimization problems

            disp('[SNLS solver: Approximate Trust Region Method]')
            disp(['Current Phase Number: ',num2str(obj.phase)])
            fprintf(' Iteration      f(x)        step        Lambda    Acceptance\n');
            formatstr = ' %5.0f        %-10.3g   %-10.3g   %-10.3g    %s';
            x0 = obj.opt.x0;
            lambda = 10; % Initial Trust Region radius
            [res,jac] = obj.cost_func(x0);
            prev_cost = res' * res;
            str = sprintf(formatstr,0,prev_cost,norm(x0),lambda,'Init');
            disp(str)

            i = 1;

            eta = obj.opt.options.LM.eta;
            Lu = obj.opt.options.LM.Lu;
            Ld = obj.opt.options.LM.Ld;
            
            while true
                A_info = jac' * jac; b = -jac' * res; 
                n = size(A_info,1);
                A_added = lambda * spdiags(diag(A_info),0,n,n);
                A = A_info + A_added;
                x0 = A \ b; % Gauss-Newton Step
                
                dummy_states = obj.states; % need to check if change in dummy state changes obj.states
                dummy_bias = obj.bias;
                dummy_imu = obj.imu;
                if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
                    dummy_arc_segments = obj.map.arc_segments;
                    dummy_assocL = obj.map.assocL;
                    dummy_assocR = obj.map.assocR;
                end

                [res_,jac_] = obj.cost_func(x0);
                cost = res_' * res_;

                ared = prev_cost - cost;
                pred = x0' * (A_added * x0 + b);
                rho = ared/pred; 
                
                if rho >= eta
                    % Current Step Accepted, update loop variables
                    res = res_;
                    jac = jac_;
                    prev_cost = cost;
                    string = 'Accepted';
                    flag = true;
                else
                    % Current Step Rejected, recover states, bias, 
                    % imu(preintegrated terms) using dummy variables
                    obj.states = dummy_states;
                    obj.bias = dummy_bias;
                    obj.imu = dummy_imu;
                    if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
                        obj.map.states = dummy_states;
                        obj.map.arc_segments = dummy_arc_segments;
                        obj.map.assocL = dummy_assocL;
                        obj.map.assocR = dummy_assocR;
                    end
                    string = 'Rejected';
                    flag = false;                    
                end
                
                step_size = norm(x0);
                
                str = sprintf(formatstr,i,cost,step_size,lambda,string);
                disp(str)
                
                lambda = obj.UpdateLambda(lambda,Lu,Ld,flag);
                
                if flag % Check Ending Criterion for Accepted Steps
                    obj.opt.flags = [];
                    obj.opt.flags = [obj.opt.flags abs(ared) > obj.opt.options.CostThres];
                    obj.opt.flags = [obj.opt.flags step_size > obj.opt.options.StepThres];
                    obj.opt.flags = [obj.opt.flags i < obj.opt.options.IterThres];
    
                    if length(find(obj.opt.flags)) ~= length(obj.opt.flags) % If any of the criterion is not met, end loop
                        
                        obj.retract(x0,'final');
    
                        disp('[Optimization Finished...]')
                        idx = find(~obj.opt.flags,1);
                        
                        if idx == 1
                            disp(['Current cost difference ',num2str(abs(ared)),' is below threshold: ',num2str(obj.opt.options.CostThres)])
                        elseif idx == 2
                            disp(['Current step size ',num2str(step_size),' is below threshold: ',num2str(obj.opt.options.StepThres)])
                        elseif idx == 3
                            disp(['Current iteration number ',num2str(i),' is above threshold: ',num2str(obj.opt.options.IterThres)])
                        end
    
                        break;
                    end

                    i = i + 1;
                end
            end

        end

        %% Trust Region Method
        function obj = TrustRegion(obj)
            % Indefinite Gauss-Newton-Powell's Dog-Leg algorithm
            % Implemented "RISE: An Incremental Trust Region Method for Robust Online Sparse Least-Squares Estimation"
            % by David M. Rosen etal
            % 
            % Implemented in MATLAB by JinHwan Jeon, 2022
            %
            % Original paper considers case for rank-deficient Jacobians,
            % but in this sensor fusion framework, Jacobian must always be
            % full rank. Therefore, 'rank-deficient Jacobian' part of the
            % algorithm is not implemented.
            %
            % Moreover, the core part of the algorithm is "incremental"
            % trust region method, but here the optimization is done in a
            % batch manner. Therefore incremental solver framework is not
            % adopted in this implementation.
            
            disp('[SNLS solver: Approximate Trust Region Method]')
            disp(['Current Phase Number: ',num2str(obj.phase)])
            fprintf(' Iteration      f(x)          step        TR_radius    Acceptance\n');
            formatstr = ' %5.0f        %-12.6g   %-10.3g   %-10.3g    %s';
            x0 = obj.opt.x0;
            tr_rad = 10; % Initial Trust Region radius
            [res,jac] = obj.cost_func(x0);
            prev_cost = res' * res;
            str = sprintf(formatstr,0,prev_cost,norm(x0),tr_rad,'Init');
            disp(str)

            i = 1;

            eta1 = obj.opt.options.TR.eta1;
            eta2 = obj.opt.options.TR.eta2;
            gamma1 = obj.opt.options.TR.gamma1;
            gamma2 = obj.opt.options.TR.gamma2;
            
            while true
                A = jac' * jac; b = -jac' * res; 
                
%                 t_start = tic;
                h_gn = A \ b; % Gauss-Newton Step
%                 elapsed_t = toc(t_start);
%                 disp(['A \ b Solver Elapsed Time: ',num2str(elapsed_t),'s'])

                alpha = (b' * b)/(b' * A * b);
                h_gd = alpha * b; % Gradient Descent Step
                
                if ~isempty(find(isnan(h_gn),1)) || ~isempty(find(isnan(h_gd),1))
                    obj.opt.info = A;
                    obj.errorAnalysis();
                end

                x0 = obj.ComputeDogLeg(h_gn,h_gd,tr_rad);
                
                dummy_states = obj.states; 
                dummy_bias = obj.bias;
                dummy_imu = obj.imu;

                if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
                    dummy_arc_segments = obj.map.arc_segments;
                    dummy_assocL = obj.map.assocL;
                    dummy_assocR = obj.map.assocR;
                end
                
                [res_,jac_] = obj.cost_func(x0);
                cost = res_' * res_;

                ared = prev_cost - cost;
                pred = b' * x0 - 1/2 * x0' * A * x0;
                rho = ared/pred; 
                
                if rho >= eta1
                    % Current Step Accepted, update loop variables
                    res = res_;
                    jac = jac_;
                    prev_cost = cost;
                    string = 'Accepted';
                    flag = true;
                else
                    % Current Step Rejected, recover states, bias, 
                    % imu(preintegrated terms) using dummy variables
                    obj.states = dummy_states;                                    
                    obj.bias = dummy_bias;
                    obj.imu = dummy_imu;
                    if strcmp(obj.mode,'2-phase') || strcmp(obj.mode,'full')
                        obj.map.states = dummy_states;
                        obj.map.arc_segments = dummy_arc_segments;
                        obj.map.assocL = dummy_assocL;
                        obj.map.assocR = dummy_assocR;
                    end
                    string = 'Rejected';
                    flag = false;
                end
                
                step_size = norm(x0);
                
                str = sprintf(formatstr,i,cost,step_size,tr_rad,string);
                disp(str)
                
                tr_rad = obj.UpdateDelta(rho,tr_rad,eta1,eta2,gamma1,gamma2);
                
                if flag % Check Ending Criterion for Accepted Steps
                    obj.opt.flags = [];
                    obj.opt.flags = [obj.opt.flags abs(ared) > obj.opt.options.CostThres];
                    obj.opt.flags = [obj.opt.flags step_size > obj.opt.options.StepThres];
                    obj.opt.flags = [obj.opt.flags i < obj.opt.options.IterThres];
                    obj.opt.flags = [obj.opt.flags tr_rad > obj.opt.options.TR.thres];

                    if length(find(obj.opt.flags)) ~= length(obj.opt.flags) % If any of the criterion is not met, end loop
                        
                        obj.retract(x0,'final');
    
                        disp('[Optimization Finished...]')
                        idx = find(~obj.opt.flags,1);
                        
                        if idx == 1
                            disp(['Current cost difference ',num2str(abs(ared)),' is below threshold: ',num2str(obj.opt.options.CostThres)])
                        elseif idx == 2
                            disp(['Current step size ',num2str(step_size),' is below threshold: ',num2str(obj.opt.options.StepThres)])
                        elseif idx == 3
                            disp(['Current iteration number ',num2str(i),' is above threshold: ',num2str(obj.opt.options.IterThres)])
                        elseif idx == 4
                            disp(['Current trust region radius ',num2str(tr_rad),' is below threshold: ',num2str(obj.opt.options.TR.thres)])
                        end
    
                        break;
                    end

                    i = i + 1;
                else
                    if tr_rad < obj.opt.options.TR.thres
                        obj.retract(x0,'final');
                        disp(['Current trust region radius ',num2str(tr_rad),' is below threshold: ',num2str(obj.opt.options.TR.thres)])
                        break;
                    end
                end
            end


            % Sparse Inverse to extract lane covariance
            obj.opt.info = jac' * jac;
            disp('[Computing Sparse Inverse of Information Matrix...]')
            obj.opt.cov = sparseinv(obj.opt.info);
            obj.extractLaneCov();
        end

        %% Optimization Cost Function
        function [res,jac] = cost_func(obj,x0)
            obj.retract(x0,'normal');
            [Pr_res,Pr_jac] = obj.CreatePrBlock();
            [MM_res,MM_jac] = obj.CreateMMBlock();
            [GNSS_res,GNSS_jac] = obj.CreateGNSSBlock();
            [WSS_res,WSS_jac] = obj.CreateWSSBlock();
            
            % Check Jacobian computation time
%             t_start = tic;
            [AS_res,AS_jac] = obj.CreateASBlock();
%             elapsed_t = toc(t_start);
%             disp(['AS Jacobian Computation Time: ',num2str(elapsed_t),'s'])
            [AS2_res,AS2_jac] = obj.CreateAS2Block();
%             [AL_res,AL_jac] = obj.CreateALBlock();
            
            res = vertcat(Pr_res,MM_res,GNSS_res,WSS_res,AS_res,AS2_res);
            jac = vertcat(Pr_jac,MM_jac,GNSS_jac,WSS_jac,AS_jac,AS2_jac);   
            obj.opt.jac = jac;
        end
        
        %% Prior Residual and Jacobian
        function [Pr_res,Pr_jac] = CreatePrBlock(obj)
            n = length(obj.states);
            blk_width = 15*n;
            blk_height = 15;
            adder = 0;

            if strcmp(obj.mode,'partial')
                blk_width = blk_width + n;
                adder = n ;                
                blk_height = blk_height + n;
            elseif strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                % Add prev_num states
                blk_height = blk_height + n;
                adder = n;
                blk_width = length(obj.opt.x0);
            end

            Veh_cov = blkdiag(obj.covs.prior.R,obj.covs.prior.V,obj.covs.prior.P);
            Bias_cov = blkdiag(obj.covs.prior.Bg,obj.covs.prior.Ba);
            
            I = zeros(1,9*3 + 6 + adder); J = I; V = I;

            % Vehicle Prior
            R_init = obj.opt.init_state.R;
            V_init = obj.opt.init_state.V;
            P_init = obj.opt.init_state.P;
            R1 = obj.states{1}.R;
            V1 = obj.states{1}.V;
            P1 = obj.states{1}.P;
            
            resR = Log_map(R_init' * R1);
            resV = R_init'* (V1 - V_init);
            resP = R_init'* (P1 - P_init);
            
            Veh_res = InvMahalanobis([resR;resV;resP],Veh_cov);

            [I_1,J_1,V_1] = sparseFormat(1:3,1:3,InvMahalanobis(InvRightJac(resR),obj.covs.prior.R));
            [I_2,J_2,V_2] = sparseFormat(4:6,4:6,InvMahalanobis(R_init',obj.covs.prior.V));
            [I_3,J_3,V_3] = sparseFormat(7:9,7:9,InvMahalanobis(R_init'*R1,obj.covs.prior.P));
            I(1:9*3) = [I_1, I_2, I_3];
            J(1:9*3) = [J_1, J_2, J_3];
            V(1:9*3) = [V_1, V_2, V_3];

            % Bias Prior
            resBg = obj.bias{1}.bg + obj.bias{1}.bgd - obj.opt.init_state.bg - obj.opt.init_state.bgd;
            resBa = obj.bias{1}.ba + obj.bias{1}.bad - obj.opt.init_state.ba - obj.opt.init_state.bad;
            
            Bias_res = InvMahalanobis([resBg;resBa],Bias_cov);
            V_1b = diag(InvMahalanobis(eye(3),obj.covs.prior.Bg));
            V_2b = diag(InvMahalanobis(eye(3),obj.covs.prior.Ba));

            I(9*3+1:9*3+6) = 9+1:9+6;
            J(9*3+1:9*3+6) = 9*n+1:9*n+6;
            V(9*3+1:9*3+6) = [V_1b V_2b];
            
            WSS_res = []; 
            if ~strcmp(obj.mode,'basic')
                % WSF Prior
                WSS_res = zeros(n,1);
                for i=1:n
                    WSF = obj.states{i}.WSF;
                    WSS_res(i) = InvMahalanobis(WSF - 1,obj.covs.prior.WSF);
                    I(9*3+6+i) = 15+i;
                    J(9*3+6+i) = 15*n+i;
                    V(9*3+6+i) = InvMahalanobis(1,obj.covs.prior.WSF);
                end       
            end            
            Pr_res = [Veh_res; Bias_res; WSS_res];
            Pr_jac = sparse(I,J,V,blk_height,blk_width);
        end
        
        %% IMU Residual and Jacobian
        function [MM_res,MM_jac] = CreateMMBlock(obj)
            n = length(obj.states);
            blk_width = 15*n;
%             blk_height = 15 * (n-1);
            grav = [0;0;-9.81];

            if ~strcmp(obj.mode,'basic')
                blk_width = blk_width + n;
                I_s = zeros(1,2*(n-1)); J_s = I_s; V_s = I_s;
                s_cov = obj.covs.imu.ScaleFactorNoise;
                WSF_res = zeros(n-1,1);

                if strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                    blk_width = length(obj.opt.x0);
                end
            end
            
            I_v = zeros(1,9*24*(n-1)); J_v = I_v; V_v = I_v;
            I_b = zeros(1,12*(n-1)); J_b = I_b; V_b = I_b;
            
            bg_cov = obj.covs.imu.GyroscopeBiasNoise;
            ba_cov = obj.covs.imu.AccelerometerBiasNoise;
            
            Veh_res = zeros(9*(n-1),1);
            Bias_res = zeros(6*(n-1),1);
            
            for i=1:n-1
                Ri = obj.states{i}.R; Rj = obj.states{i+1}.R;
                Vi = obj.states{i}.V; Vj = obj.states{i+1}.V;
                Pi = obj.states{i}.P; Pj = obj.states{i+1}.P;
                
                bgi = obj.bias{i}.bg; bgj = obj.bias{i+1}.bg;
                bai = obj.bias{i}.ba; baj = obj.bias{i+1}.ba;
                bgdi = obj.bias{i}.bgd; bgdj = obj.bias{i+1}.bgd;
                badi = obj.bias{i}.bad; badj = obj.bias{i+1}.bad;

                % Use pre-integrated values
                JdelRij_bg = obj.imu{i}.JdelRij_bg;
                JdelVij_bg = obj.imu{i}.JdelVij_bg;
                JdelVij_ba = obj.imu{i}.JdelVij_ba;
                JdelPij_bg = obj.imu{i}.JdelPij_bg;
                JdelPij_ba = obj.imu{i}.JdelPij_ba;

                delRij = obj.imu{i}.delRij;
                delVij = obj.imu{i}.delVij;
                delPij = obj.imu{i}.delPij;
                
                dtij = obj.imu{i}.dtij;

                Covij = obj.imu{i}.Covij;
                
                % IMU residual and jacobian
                resR = Log_map((delRij * Exp_map(JdelRij_bg * bgdi))' * Ri' * Rj);
                resV = Ri' * (Vj - Vi - grav *dtij) - (delVij + JdelVij_bg * bgdi + JdelVij_ba * badi);
                resP = Ri' * (Pj - Pi - Vi * dtij - 1/2 * grav * dtij^2) - (delPij + JdelPij_bg * bgdi + JdelPij_ba * badi);

                Veh_res((9*(i-1)+1):9*i) = InvMahalanobis([resR;resV;resP],Covij);

                [Ji,Jj,Jb] = obj.getIMUJac(resR,Ri,Rj,Vi,Vj,Pi,Pj,bgdi,JdelRij_bg,...
                                           JdelVij_bg,JdelVij_ba,JdelPij_bg,JdelPij_ba,dtij);
                
                [I_1,J_1,V_1] = sparseFormat(9*(i-1)+1:9*(i-1)+9,9*(i-1)+1:9*(i-1)+9,InvMahalanobis(Ji,Covij));
                [I_2,J_2,V_2] = sparseFormat(9*(i-1)+1:9*(i-1)+9,9*i+1:9*i+9,InvMahalanobis(Jj,Covij));
                [I_3,J_3,V_3] = sparseFormat(9*(i-1)+1:9*(i-1)+9,9*n+6*(i-1)+1:9*n+6*(i-1)+6,InvMahalanobis(Jb,Covij));

                I_v((9*24*(i-1)+1):9*24*i) = [I_1,I_2,I_3];
                J_v((9*24*(i-1)+1):9*24*i) = [J_1,J_2,J_3];
                V_v((9*24*(i-1)+1):9*24*i) = [V_1,V_2,V_3];
                
                % Bias residual and jacobian
                resBg = bgj + bgdj - (bgi + bgdi);
                resBa = baj + badj - (bai + badi);

                Bias_res((6*(i-1)+1):6*i) = InvMahalanobis([resBg;resBa],blkdiag(dtij * bg_cov, dtij * ba_cov));

                I_b1 = 6*(i-1)+1:6*(i-1)+3; J_b1 = 9*n+6*(i-1)+1:9*n+6*(i-1)+3; V_b1 = diag(InvMahalanobis(-eye(3),dtij * bg_cov));
                I_b2 = 6*(i-1)+1:6*(i-1)+3; J_b2 = 9*n+6*i+1:9*n+6*i+3; V_b2 = diag(InvMahalanobis(eye(3),dtij * bg_cov));
                I_b3 = 6*(i-1)+4:6*(i-1)+6; J_b3 = 9*n+6*(i-1)+4:9*n+6*(i-1)+6; V_b3 = diag(InvMahalanobis(-eye(3),dtij * ba_cov));
                I_b4 = 6*(i-1)+4:6*(i-1)+6; J_b4 = 9*n+6*i+4:9*n+6*i+6; V_b4 = diag(InvMahalanobis(eye(3),dtij * ba_cov));

                I_b((12*(i-1)+1):12*i) = [I_b1 I_b2 I_b3 I_b4];
                J_b((12*(i-1)+1):12*i) = [J_b1 J_b2 J_b3 J_b4];
                V_b((12*(i-1)+1):12*i) = [V_b1 V_b2 V_b3 V_b4];
                
                if ~strcmp(obj.mode,'basic')
                    % WSF residual and jacobian
                    Si = obj.states{i}.WSF; Sj = obj.states{i+1}.WSF;
                    resS = Sj - Si;
                    WSF_res(i) = InvMahalanobis(resS,dtij * s_cov);
                    I_s(2*i-1:2*i) = [i,i];
                    J_s(2*i-1:2*i) = 15*n+i:15*n+i+1;
                    V_s(2*i-1:2*i) = InvMahalanobis([-1,1],dtij * s_cov);
                end
            end
            Veh_jac = sparse(I_v,J_v,V_v,9*(n-1),blk_width);
            Bias_jac = sparse(I_b,J_b,V_b,6*(n-1),blk_width);

            MM_res = [Veh_res;Bias_res];
            MM_jac = [Veh_jac;Bias_jac];
            
            if ~strcmp(obj.mode,'basic')
                % Augment WSF residual and jacobian
                WSF_jac = sparse(I_s,J_s,V_s,n-1,blk_width);
                
                MM_res = [MM_res; WSF_res];
                MM_jac = [MM_jac; WSF_jac];
            end

        end

        %% GNSS Residual and Jacobian
        function [GNSS_res,GNSS_jac] = CreateGNSSBlock(obj)
            idxs = obj.gnss.state_idxs;
            n = length(idxs);
            
            m = length(obj.states);
            blk_width = 15*m;

            if strcmp(obj.mode,'partial')
                blk_width = blk_width + m;
            elseif strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                blk_width = length(obj.opt.x0);
            end
            blk_height = 3*n;
            GNSS_res = zeros(blk_height,1);
            I_g = zeros(1,3*3*n); J_g = I_g; V_g = I_g;

            for i=1:n
                idx = idxs(i); % state_idx
                Ri = obj.states{idx}.R;
                Pi = obj.states{idx}.P;
                P_meas = lla2enu(obj.gnss.pos(i,:),obj.gnss.lla0,'ellipsoid')';
                hAcc = obj.gnss.hAcc(i);
                vAcc = obj.gnss.vAcc(i);
                cov = diag([hAcc^2 hAcc^2 vAcc^2]);

                GNSS_res((3*(i-1)+1):3*i) = InvMahalanobis(Pi-P_meas,cov);
                
                [I_,J_,V_] = sparseFormat((3*(i-1)+1):3*i,9*(idx-1)+7:9*(idx-1)+9,InvMahalanobis(Ri,cov));
                I_g((9*(i-1)+1):9*i) = I_;
                J_g((9*(i-1)+1):9*i) = J_;
                V_g((9*(i-1)+1):9*i) = V_;
            end
            GNSS_jac = sparse(I_g,J_g,V_g,blk_height,blk_width);
        end
        
        %% WSS Residual and Jacobian
        function [WSS_res,WSS_jac] = CreateWSSBlock(obj)
            if strcmp(obj.mode,'basic')
                WSS_res = [];
                WSS_jac = [];
            else
                n = length(obj.states);
                blk_width = 16*n;
                
                if strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                    blk_width = length(obj.opt.x0);
                end
                
                L = obj.params;
    
                WSS_res = zeros(3*(n-2),1);
                m = 9*3+3;
                I = zeros(1,m*(n-2)); J = I; V = I; 
                wss_cov = obj.covs.wss;

                for i=2:n-1
                    Ri = obj.states{i}.R; Vi = obj.states{i}.V; 
                    Si = obj.states{i}.WSF;
                    wi = obj.imu{i}.gyro(1);
                    
                    bgi = obj.bias{i}.bg; bgdi = obj.bias{i}.bgd; 
                    
                    idx = obj.can.state_idxs(i);
                    V_meas = [obj.can.whl_spd(idx); 0; 0]; % Vehicle Nonholonomic Constraint
                    
                    res =  1/Si * (Ri' * Vi - skew(wi - bgi - bgdi) * L) - V_meas;
                    WSS_res(3*(i-2)+1:3*(i-2)+3) = InvMahalanobis(res,wss_cov);
    
                    [Jrv,Jbg,Js] = obj.getWSSJac(Ri,Vi,Si,wi,bgi,bgdi,L);

                    [I_rv,J_rv,V_rv] = sparseFormat(3*(i-2)+1:3*(i-2)+3,9*(i-1)+1:9*(i-1)+6,InvMahalanobis(Jrv,wss_cov));
                    [I_bg,J_bg,V_bg] = sparseFormat(3*(i-2)+1:3*(i-2)+3,9*n+6*(i-1)+1:9*n+6*(i-1)+3,InvMahalanobis(Jbg,wss_cov));
                    [I_s,J_s,V_s] = sparseFormat(3*(i-2)+1:3*(i-2)+3,15*n+i,InvMahalanobis(Js,wss_cov));
                    
                    I(m*(i-2)+1:m*(i-1)) = [I_rv I_bg I_s];
                    J(m*(i-2)+1:m*(i-1)) = [J_rv J_bg J_s];
                    V(m*(i-2)+1:m*(i-1)) = [V_rv V_bg V_s];
                end

                WSS_jac = sparse(I,J,V,3*(n-2),blk_width);
            end
        end
        
        %% Arc Spline based Measurement Residual and Jacobian 
        function [AS_res,AS_jac] = CreateASBlock(obj)
            % Arc Spline based lane measurement model, Phase 2
            % Pre-computed data association matrix is used for state-map
            % segment matching
            if ~strcmp(obj.mode,'2-phase') && ~strcmp(obj.mode,'full')
                AS_res = []; AS_jac = [];
            else
                % Due to the complexity of computing jacobian of
                % measurement function algebraically, numerical jacobian
                % computation is implemented for this block
                % For computation speed, forward-difference is adopted
                % 
                % This part is to be merged with phase 3 functions after
                % debugging process

                n = size(obj.lane.FactorValidIntvs,1);
                m = length(obj.states);
                
                blk_heightL = nnz(obj.map.assocL);
                blk_heightR = nnz(obj.map.assocR);
                blk_width = length(obj.opt.x0);
                AS_resL = zeros(blk_heightL,1);
                AS_resR = zeros(blk_heightR,1);
                
                I_L = zeros(1,9*blk_heightL + 2 * sum(obj.map.assocL,'all'));
                J_L = I_L; V_L = I_L;

                I_R = zeros(1,9*blk_heightR + 2 * sum(obj.map.assocR,'all'));
                J_R = I_R; V_R = I_R;

                cntL = 1; cntR = 1;
                jac_cntL = 0; jac_cntR = 0;
                for i=1:n
                    lb = obj.lane.FactorValidIntvs(i,1);
                    ub = obj.lane.FactorValidIntvs(i,2);

                    leftSegIdx = obj.map.segment_info(1,i);
                    rightSegIdx = obj.map.segment_info(2,i);
                    leftSeg = obj.map.arc_segments{leftSegIdx};
                    rightSeg = obj.map.arc_segments{rightSegIdx};

                    for j=lb:ub
                        R = obj.states{j}.R;
                        P = obj.states{j}.P;
                        lane_idx = obj.lane.state_idxs(j);
                        
%                         if obj.phase == 1
%                             max_prev = 1;
%                         elseif obj.phase == 2
%                             max_prev = ;
%                         end

                        for k=1:obj.lane.prev_num                           
                            leftSubSegIdx = obj.map.assocL(j,k);
                            rightSubSegIdx = obj.map.assocR(j,k);
    
                            if leftSubSegIdx > 0
                                initParams = [leftSeg.x0, leftSeg.y0, leftSeg.tau0];
                                [r_jac, p_jac, param_jac, kappa_jac, L_jac, anchored_res] = ...
                                obj.getASJac(R,P,initParams,leftSeg.kappa,leftSeg.L,leftSubSegIdx,lane_idx,k,'left');
    
                                AS_resL(cntL) = anchored_res;
    
                                % Warning sign for complex residual
                                if ~isreal(anchored_res)
                                    warning(['Complex residual(Left) detected at state_idx: ',num2str(j),...
                                             ', segment idx: ',num2str(leftSegIdx),...
                                             ', sub-segment idx: ',num2str(leftSubSegIdx)])                                
                                end
    
                                idxL = 16*m+obj.opt.seg_tracker(leftSegIdx);
                                idxkL = idxL + 3;
                                idxlL = idxL + 3 + length(leftSeg.kappa);
                                [I_rL,J_rL,V_rL] = sparseFormat(cntL,9*(j-1)+1:9*(j-1)+3,r_jac);
                                [I_pL,J_pL,V_pL] = sparseFormat(cntL,9*(j-1)+7:9*(j-1)+9,p_jac);
                                [I_PL,J_PL,V_PL] = sparseFormat(cntL,idxL+1:idxL+3,param_jac);
                                [I_kL,J_kL,V_kL] = sparseFormat(cntL,idxkL+1:idxkL+leftSubSegIdx,kappa_jac);
                                [I_lL,J_lL,V_lL] = sparseFormat(cntL,idxlL+1:idxlL+leftSubSegIdx,L_jac);
                                
                                I_L(jac_cntL+1:jac_cntL+9+2*leftSubSegIdx) = [I_rL, I_pL, I_PL, I_kL, I_lL];
                                J_L(jac_cntL+1:jac_cntL+9+2*leftSubSegIdx) = [J_rL, J_pL, J_PL, J_kL, J_lL];
                                V_L(jac_cntL+1:jac_cntL+9+2*leftSubSegIdx) = [V_rL, V_pL, V_PL, V_kL, V_lL];
    
                                cntL = cntL + 1;
                                jac_cntL = jac_cntL + 9 + 2 * leftSubSegIdx;
                            end
    
                            if rightSubSegIdx > 0
                                initParams = [rightSeg.x0, rightSeg.y0, rightSeg.tau0];
                                [r_jac, p_jac, param_jac, kappa_jac, L_jac, anchored_res] = ...
                                obj.getASJac(R,P,initParams,rightSeg.kappa,rightSeg.L,rightSubSegIdx,lane_idx,k,'right');
    
                                AS_resR(cntR) = anchored_res;
    
                                % Warning sign for complex residual
                                if ~isreal(anchored_res)
                                    warning(['Complex residual(Right) detected at state_idx: ',num2str(j),...
                                             ', segment idx: ',num2str(rightSegIdx),...
                                             ', sub-segment idx: ',num2str(rightSubSegIdx)])                                
                                end
    
                                idxR = 16*m+obj.opt.seg_tracker(rightSegIdx);
                                idxkR = idxR + 3;
                                idxlR = idxR + 3 + length(rightSeg.kappa);
                                [I_rR,J_rR,V_rR] = sparseFormat(cntR,9*(j-1)+1:9*(j-1)+3,r_jac);
                                [I_pR,J_pR,V_pR] = sparseFormat(cntR,9*(j-1)+7:9*(j-1)+9,p_jac);
                                [I_PR,J_PR,V_PR] = sparseFormat(cntR,idxR+1:idxR+3,param_jac);
                                [I_kR,J_kR,V_kR] = sparseFormat(cntR,idxkR+1:idxkR+rightSubSegIdx,kappa_jac);
                                [I_lR,J_lR,V_lR] = sparseFormat(cntR,idxlR+1:idxlR+rightSubSegIdx,L_jac);
                                
                                I_R(jac_cntR+1:jac_cntR+9+2*rightSubSegIdx) = [I_rR, I_pR, I_PR, I_kR, I_lR];
                                J_R(jac_cntR+1:jac_cntR+9+2*rightSubSegIdx) = [J_rR, J_pR, J_PR, J_kR, J_lR];
                                V_R(jac_cntR+1:jac_cntR+9+2*rightSubSegIdx) = [V_rR, V_pR, V_PR, V_kR, V_lR];
    
                                cntR = cntR + 1;
                                jac_cntR = jac_cntR + 9 + 2 * rightSubSegIdx;
                            end
                        end
                    end
                end
                AS_res = vertcat(AS_resL,AS_resR);
                AS_jacL = sparse(I_L,J_L,V_L,blk_heightL,blk_width);
                AS_jacR = sparse(I_R,J_R,V_R,blk_heightR,blk_width);
                AS_jac = vertcat(AS_jacL,AS_jacR);

                obj.opt.AS_res = AS_res;
                obj.opt.AS_jac = AS_jac;
            end
        end

        %% Arc Spline based Measurement 2 Residual and Jacobian Phase 2
        function [AS2_res,AS2_jac] = CreateAS2Block(obj)
            if ~strcmp(obj.mode,'2-phase') && ~strcmp(obj.mode,'full')
                AS2_res = []; AS2_jac = [];
            else
                n = length(obj.map.arc_segments);
                m = length(obj.states);
%                 blk_height = 2*n;
%                 blk_width = length(obj.opt.x0);
%                 AS2_resF = zeros(blk_height,1);
%                 AS2_resB = zeros(blk_height,1);
%                 
%                 I_F = []; J_F = []; V_F = [];
%                 I_B = []; J_B = []; V_B = [];
% 
%                 % First point anchoring equation
%                 for i=1:n
%                     [intvIdx,dir] = obj.minColIdx(obj.map.segment_info,i);
%                     state_idx = obj.lane.FactorValidIntvs(intvIdx,1);
%                     R = obj.states{state_idx}.R; P = obj.states{state_idx}.P;
%                     lane_idx = obj.lane.state_idxs(state_idx);
%                     seg = obj.map.arc_segments{i};
%                     initParams = [seg.x0, seg.y0, seg.tau0];
%                     kappa = seg.kappa; L = seg.L;
%                     [r_jac,p_jac,param_jac,kappa_jac,L_jac,anchored_res] = ...
%                     obj.getAS2Jac(R,P,initParams,kappa,L,lane_idx,dir,'init');
% 
%                     AS2_resF(2*i-1:2*i) = anchored_res;
% 
%                     idxF = 16 * m + obj.opt.seg_tracker(i);
%                     idxkF = idxF + 3;
%                     idxlF = idxF + 3 + length(seg.kappa);
%                     [I_rF,J_rF,V_rF] = sparseFormat(2*i-1:2*i,9*(state_idx-1)+1:9*(state_idx-1)+3,r_jac);
%                     [I_pF,J_pF,V_pF] = sparseFormat(2*i-1:2*i,9*(state_idx-1)+7:9*(state_idx-1)+9,p_jac);
%                     [I_PF,J_PF,V_PF] = sparseFormat(2*i-1:2*i,idxF+1:idxF+3,param_jac);
%                     % Jacobian for subsegment parameters should
%                     % be 0 algebraically
%                     [I_kF,J_kF,V_kF] = sparseFormat(2*i-1:2*i,idxkF+1:idxkF+length(seg.kappa),kappa_jac);
%                     [I_lF,J_lF,V_lF] = sparseFormat(2*i-1:2*i,idxlF+1:idxlF+length(seg.L),L_jac);
%                     I_F = [I_F I_rF I_pF I_PF I_kF I_lF];
%                     J_F = [J_F J_rF J_pF J_PF J_kF J_lF];
%                     V_F = [V_F V_rF V_pF V_PF V_kF V_lF];
%                 end
% 
%                 % Final point anchoring equation
%                 for i=1:n
%                     [intvIdx,dir] = obj.maxColIdx(obj.map.segment_info,i);
%                     state_idx = obj.lane.FactorValidIntvs(intvIdx,2);
%                     R = obj.states{state_idx}.R; P = obj.states{state_idx}.P;
%                     lane_idx = obj.lane.state_idxs(state_idx);
%                     seg = obj.map.arc_segments{i};
%                     initParams = [seg.x0, seg.y0, seg.tau0];
%                     kappa = seg.kappa; L = seg.L;
%                     [r_jac,p_jac,param_jac,kappa_jac,L_jac,anchored_res] = ...
%                     obj.getAS2Jac(R,P,initParams,kappa,L,lane_idx,dir,'last');
% 
%                     AS2_resB(2*i-1:2*i) = anchored_res;
%                     
%                     idxB = 16 * m + obj.opt.seg_tracker(i);
%                     idxkB = idxB + 3;
%                     idxlB = idxB + 3 + length(kappa);
%                     [I_rB,J_rB,V_rB] = sparseFormat(2*i-1:2*i,9*(state_idx-1)+1:9*(state_idx-1)+3,r_jac);
%                     [I_pB,J_pB,V_pB] = sparseFormat(2*i-1:2*i,9*(state_idx-1)+7:9*(state_idx-1)+9,p_jac);
%                     [I_PB,J_PB,V_PB] = sparseFormat(2*i-1:2*i,idxB+1:idxB+3,param_jac);
%                     [I_kB,J_kB,V_kB] = sparseFormat(2*i-1:2*i,idxkB+1:idxkB+length(kappa),kappa_jac);
%                     [I_lB,J_lB,V_lB] = sparseFormat(2*i-1:2*i,idxlB+1:idxlB+length(L),L_jac);
%                     I_B = [I_B I_rB I_pB I_PB I_kB I_lB];
%                     J_B = [J_B J_rB J_pB J_PB J_kB J_lB];
%                     V_B = [V_B V_rB V_pB V_PB V_kB V_lB];
%                 end
% 
%                 AS2_res = vertcat(AS2_resF,AS2_resB);
%                 AS2_jacF = sparse(I_F,J_F,V_F,blk_height,blk_width);
%                 AS2_jacB = sparse(I_B,J_B,V_B,blk_height,blk_width);
%                 AS2_jac = vertcat(AS2_jacF,AS2_jacB);
% 
%                 obj.opt.AS2_jac = AS2_jac;
%                 obj.opt.AS2_res = AS2_res;

                numSubSeg = 0;
                I_size = 0;
                for i=1:n
                    % # of Node Points = # of Sub-segments + 1
                    N = length(obj.map.arc_segments{i}.kappa);
                    numSubSeg = numSubSeg + N + 1;
                    I_size = I_size + 2*N^2 + 20 * N + 18;
                end
                blk_height = 2 * numSubSeg;
                blk_width = length(obj.opt.x0);
                AS2_res = zeros(blk_height,1);
                I = zeros(1,I_size); J = I; V = I;
                
                cnt = 0;
                jac_cnt = 0;
                for i=1:n
                    seg = obj.map.arc_segments{i};
                    bnds = seg.bnds;
                    [r,c] = find(obj.map.segment_info == i);
                    stateIntvs = obj.lane.FactorValidIntvs(c(1):c(end),1:2);

                    % Init Node Point
                    [stateIdx, dir] = obj.map.findStateIdx(stateIntvs,bnds(1,1),r);
                    laneIdx = obj.lane.state_idxs(stateIdx);
                    R = obj.states{stateIdx}.R;
                    P = obj.states{stateIdx}.P;
                    initParams = [seg.x0, seg.y0, seg.tau0];
                    kappa = seg.kappa; L = seg.L;

                    idxInit = 16 * m + obj.opt.seg_tracker(i);
                    idxKappa = idxInit + 3;
                    idxL = idxKappa + length(seg.kappa);

                    [r_jac,p_jac,param_jac,~,~,anchored_res] = obj.getAS2Jac(R,P,initParams,kappa,L,laneIdx,dir,0);
                    
                    AS2_res(2*cnt+1:2*cnt+2) = anchored_res;
                    
                    [I_r,J_r,V_r] = sparseFormat(2*cnt+1:2*cnt+2,9*(stateIdx-1)+1:9*(stateIdx-1)+3,r_jac);
                    [I_p,J_p,V_p] = sparseFormat(2*cnt+1:2*cnt+2,9*(stateIdx-1)+7:9*(stateIdx-1)+9,p_jac);
                    [I_P,J_P,V_P] = sparseFormat(2*cnt+1:2*cnt+2,idxInit+1:idxInit+3,param_jac);
                    
                    I(jac_cnt+1:jac_cnt+9*2) = [I_r,I_p,I_P];
                    J(jac_cnt+1:jac_cnt+9*2) = [J_r,J_p,J_P];
                    V(jac_cnt+1:jac_cnt+9*2) = [V_r,V_p,V_P];

                    cnt = cnt + 1;
                    jac_cnt = jac_cnt + 9*2;

                    % Remaining Node Points
                    for j=1:size(bnds,1)
                        [stateIdx, dir] = obj.map.findStateIdx(stateIntvs,bnds(j,2),r);
                        laneIdx = obj.lane.state_idxs(stateIdx);
                        R = obj.states{stateIdx}.R;
                        P = obj.states{stateIdx}.P;
                        initParams = [seg.x0, seg.y0, seg.tau0];
                        kappa = seg.kappa; L = seg.L;
                        [r_jac,p_jac,param_jac,kappa_jac,L_jac,anchored_res] = obj.getAS2Jac(R,P,initParams,kappa,L,laneIdx,dir,j);
                        
                        AS2_res(2*cnt+1:2*cnt+2) = anchored_res;
                        
                        [I_r,J_r,V_r] = sparseFormat(2*cnt+1:2*cnt+2,9*(stateIdx-1)+1:9*(stateIdx-1)+3,r_jac);
                        [I_p,J_p,V_p] = sparseFormat(2*cnt+1:2*cnt+2,9*(stateIdx-1)+7:9*(stateIdx-1)+9,p_jac);
                        [I_P,J_P,V_P] = sparseFormat(2*cnt+1:2*cnt+2,idxInit+1:idxInit+3,param_jac);
                        [I_k,J_k,V_k] = sparseFormat(2*cnt+1:2*cnt+2,idxKappa+1:idxKappa+j,kappa_jac);
                        [I_L,J_L,V_L] = sparseFormat(2*cnt+1:2*cnt+2,idxL+1:idxL+j,L_jac);

                        I(jac_cnt+1:jac_cnt+9*2+2*2*j) = [I_r,I_p,I_P,I_k,I_L];
                        J(jac_cnt+1:jac_cnt+9*2+2*2*j) = [J_r,J_p,J_P,J_k,J_L];
                        V(jac_cnt+1:jac_cnt+9*2+2*2*j) = [V_r,V_p,V_P,V_k,V_L];

                        cnt = cnt + 1;
                        jac_cnt = jac_cnt + 9*2 + 2*2*j;
                    end
                end

                AS2_jac = sparse(I,J,V,blk_height,blk_width);
            end

        end

        %% Retraction
        function obj = retract(obj,delta,mode)
            
            n = length(obj.states);
            state_delta = delta(1:9*n);
            bias_ddelta = delta(9*n+1:15*n); % delta of delta value!

            obj.retractStates(state_delta);
            obj.retractBias(bias_ddelta,mode);
            
            if ~strcmp(obj.mode,'basic')
                wsf_delta = delta(15*n+1:16*n);
                obj.retractWSF(wsf_delta);
                
                if strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                    arc_delta = delta(16*n+1:end); % Change this part if more models are added
                    obj.retractLane(arc_delta);
                end
            end
            
        end

        %% State Retraction
        function obj = retractStates(obj,state_delta)            
            n = length(obj.states);
            for i=1:n
                deltaR = state_delta(9*(i-1)+1:9*(i-1)+3);
                deltaV = state_delta(9*(i-1)+4:9*(i-1)+6);
                deltaP = state_delta(9*(i-1)+7:9*(i-1)+9);  

                R = obj.states{i}.R; V = obj.states{i}.V; P = obj.states{i}.P; 
                P_ = P + R * deltaP;
                V_ = V + deltaV;
                R_ = R * Exp_map(deltaR);

                obj.states{i}.R = R_;
                obj.states{i}.V = V_;
                obj.states{i}.P = P_;
            end            
        end
        
        %% Bias Retraction
        function obj = retractBias(obj,bias_ddelta,mode)            
            n = length(obj.bias);
            idxs = [];
            for i=1:n                
                bgd_ = bias_ddelta(6*(i-1)+1:6*(i-1)+3);
                bad_ = bias_ddelta(6*(i-1)+4:6*(i-1)+6);

                new_bgd = bgd_ + obj.bias{i}.bgd;
                new_bad = bad_ + obj.bias{i}.bad;
                
                % Compare with threshold and update if needed
                flag = false;
            
                if strcmp(mode,'final') 
                   % After finishing optimization, all delta terms are
                   % transferred to the original bias variables
                   % flush delta values!
                   bg_ = obj.bias{i}.bg + new_bgd;
                   bgd_ = zeros(3,1);
                   ba_ = obj.bias{i}.ba + new_bad;
                   bad_ = zeros(3,1);

                elseif strcmp(mode,'normal')
                    % If delta bias values are larger than threshold,
                    % perform repropagation (IMU)
                    % Else, accumulate delta value
                    if norm(new_bgd) > obj.bgd_thres
                        bg_ = obj.bias{i}.bg + new_bgd;
                        bgd_ = zeros(3,1);
                        flag = true;
                    else
                        bg_= obj.bias{i}.bg;
                        bgd_ = new_bgd;
                    end
    
                    if norm(new_bad) > obj.bad_thres
                        ba_ = obj.bias{i}.ba + new_bad;
                        bad_ = zeros(3,1);
                        flag = true;
                    else
                        ba_ = obj.bias{i}.ba;
                        bad_ = new_bad;
                    end
                end
                

                if flag && i~=n
                    idxs = [idxs i];
                end
                
                obj.bias{i}.bg = bg_;
                obj.bias{i}.bgd = bgd_;
                obj.bias{i}.ba = ba_;
                obj.bias{i}.bad = bad_;
            end

            if ~isempty(idxs)
                % If repropagation is needed                
                obj.integrate(idxs);
            end
        end
        
        %% WSF Retraction
        function obj = retractWSF(obj,wsf_delta)
            n = length(obj.states);
            for i=1:n
                deltaWSF = wsf_delta(i);                
                WSF = obj.states{i}.WSF;

                WSF_ = WSF + deltaWSF;                
                obj.states{i}.WSF = WSF_;
            end            
        end
        
        %% Lane Retraction
        function obj = retractLane(obj,arc_delta)
            % Classify delta values according to segment number
            cnt = 0;
            arc_delta_params = {};
            
            for i=1:length(obj.map.subseg_cnt)
                delta_params = struct();
                num_subseg = obj.map.subseg_cnt(i);
                subseg_vars = arc_delta(cnt+1:cnt+3+2*num_subseg);

                delta_params.x0 = subseg_vars(1);
                delta_params.y0 = subseg_vars(2);
                delta_params.tau0 = subseg_vars(3);
                
                kappaL = subseg_vars(4:end)';
                delta_params.kappa = kappaL(1:num_subseg);
                delta_params.L = kappaL(num_subseg+1:end);
                
                arc_delta_params = [arc_delta_params {delta_params}];
                cnt = cnt + 3 + 2 * num_subseg;
            end
            obj.map.update(obj.states,arc_delta_params);
        end
        
        %% Arc Spline based Numerical Jacobian Computation Phase 2
        function [r_jac,p_jac,param_jac,kappa_jac,L_jac,anchored_res] = getASJac(obj,R,P,initParams,kappa,L,SubSegIdx,lane_idx,k,dir)
            eps = 1e-6; eps_kappa = 1e-10; minL = obj.lane.minL;
            anchored_res = obj.getASRes(R,P,initParams,kappa,L,SubSegIdx,lane_idx,k,dir);
            
            r_jac = zeros(1,3); p_jac = zeros(1,3); 
            param_jac = zeros(1,3); kappa_jac = zeros(1,SubSegIdx);
            L_jac = zeros(1,SubSegIdx);
            if strcmp(dir,'left')
                cov = obj.lane.lystd(lane_idx,k)^2;
            elseif strcmp(dir,'right')
                cov = obj.lane.rystd(lane_idx,k)^2;
            end
            
            % Add perturbation to all variables of
            % interest to find jacobian

            % 1: Perturbation to Rotational Matrix

            for i=1:3
                rot_ptb_vec = zeros(3,1);
                rot_ptb_vec(i) = eps;
                R_ptb = R * Exp_map(rot_ptb_vec); 
                res_ptb = obj.getASRes(R_ptb,P,initParams,kappa,L,SubSegIdx,lane_idx,k,dir);
                r_jac(i) = (res_ptb - anchored_res)/eps;
            end
            r_jac = InvMahalanobis(r_jac,cov);

            % 2: Perturbation to Position
            
            for i=1:3
                pos_ptb_vec = zeros(3,1);
                pos_ptb_vec(i) = eps;
                P_ptb = P + R * pos_ptb_vec;
                res_ptb = obj.getASRes(R,P_ptb,initParams,kappa,L,SubSegIdx,lane_idx,k,dir);
                p_jac(i) = (res_ptb - anchored_res)/eps;
            end
            p_jac = InvMahalanobis(p_jac,cov);

            % 3: Perturbation to Segment initParams
            for i=1:3
                initParams_ptb_vec = zeros(1,3);
                initParams_ptb_vec(i) = eps;
                initParams_ptb = initParams + initParams_ptb_vec;
                res_ptb = obj.getASRes(R,P,initParams_ptb,kappa,L,SubSegIdx,lane_idx,k,dir);
                param_jac(i) = (res_ptb - anchored_res)/eps;
            end
            param_jac = InvMahalanobis(param_jac,cov);

            % 4: Perturbation to Sub-Segment Params (curvature)
            for i=1:SubSegIdx
                % For sub-segment parameter perturbation, variables "after"
                % the target sub-segment index are not used!
                kappa_ptb_vec = zeros(1,length(kappa));
                kappa_ptb_vec(i) = eps_kappa;
                kappa_ptb = kappa + kappa_ptb_vec;
                res_ptb = obj.getASRes(R,P,initParams,kappa_ptb,L,SubSegIdx,lane_idx,k,dir);
                kappa_jac(i) = (res_ptb - anchored_res)/eps_kappa;
            end
            kappa_jac = InvMahalanobis(kappa_jac,cov);

            % 5: Perturbation to Sub-Segment Params (arc-length)
            for i=1:SubSegIdx
                L_ptb_vec = zeros(1,length(L));
                L_ptb_vec(i) = eps;
                % L positive constraint
                % actual optimization variable is sqrt(L)
                L_ptb = ((L - minL).^(1/2) + L_ptb_vec).^2 + minL;
                res_ptb = obj.getASRes(R,P,initParams,kappa,L_ptb,SubSegIdx,lane_idx,k,dir);
                L_jac(i) = (res_ptb - anchored_res)/eps;
            end
            L_jac = InvMahalanobis(L_jac,cov);
            
            anchored_res = InvMahalanobis(anchored_res,cov);
        end
        
        %% Arc Spline based measurement residual Phase 2
        function res = getASRes(obj,R,P,initParams,kappa,L,SubSegIdx,lane_idx,k,dir)
            rpy = dcm2rpy(R); psi = rpy(3);
            R2d = [cos(psi) -sin(psi);
                   sin(psi) cos(psi)];
            if strcmp(dir,'left')
                dij = obj.lane.ly(lane_idx,k);
%                 zij = obj.lane.lz(lane_idx,k);
            elseif strcmp(dir,'right')
                dij = obj.lane.ry(lane_idx,k);
%                 zij = obj.lane.rz(lane_idx,k);
            end
            % Propgate arc center point until the matched sub-segment 
            centerPoints = obj.propCenter(initParams,kappa,L);
            Xc = centerPoints(:,SubSegIdx);
%             Pl = P + R * [10*(k-1);dij;zij];
%             
%             kappa_ = kappa(SubSegIdx);
%             res = sqrt((Xc - Pl(1:2))' * (Xc - Pl(1:2))) - abs(1/kappa_);            

            % WARNING! Need to double check if this logic is appropriate
            % Computing "distance" residual
            % Very unstable, creating complex residual often
            Pp = P(1:2) + R2d * [10*(k-1);0];
            Xc_b = R2d' * (Xc - Pp); xc_b = Xc_b(1); yc_b = Xc_b(2);
            kappa_ = kappa(SubSegIdx);
            if kappa_ > 0
                res = yc_b - sqrt((1/kappa_)^2 - xc_b^2) - dij;
            else
                res = yc_b + sqrt((1/kappa_)^2 - xc_b^2) - dij;
            end
        end

        %% Arc Spline based Numerical Jacobian Computation 2 Phase 2
        function [r_jac,p_jac,param_jac,kappa_jac,L_jac,anchored_res] = getAS2Jac(obj,R,P,initParams,kappa,L,lane_idx,dir,SubSegIdx)
            eps = 1e-6; eps_kappa = 1e-10; minL = obj.lane.minL;
%             eps_L = 1e-6 * linspace(1,length(L));
            eps_L = 5;
            anchored_res = obj.getAS2Res(R,P,initParams,kappa,L,lane_idx,dir,SubSegIdx);
            
            r_jac = zeros(2,3); p_jac = zeros(2,3); 
            param_jac = zeros(2,3);
            kappa_jac = zeros(2,SubSegIdx); L_jac = zeros(2,SubSegIdx);

            if dir == 1
                cov = diag([0.01^2, obj.lane.lystd(lane_idx,1)^2]);
            elseif dir == 2
                cov = diag([0.01^2, obj.lane.rystd(lane_idx,1)^2]);
            end
%             cov = 1e-5 * eye(1);
            
            % Add perturbation to all variables of
            % interest to find jacobian

            % 1: Perturbation to Rotational Matrix
            
            for i=1:3
                rot_ptb_vec = zeros(3,1);
                rot_ptb_vec(i) = eps;
                R_ptb = R * Exp_map(rot_ptb_vec);
                res_ptb = obj.getAS2Res(R_ptb,P,initParams,kappa,L,lane_idx,dir,SubSegIdx);               
                r_jac(:,i) = (res_ptb - anchored_res)/eps;
            end
            r_jac = InvMahalanobis(r_jac,cov);

            % 2: Perturbation to Position
            
            for i=1:3
                pos_ptb_vec = zeros(3,1);
                pos_ptb_vec(i) = eps;
                P_ptb = P + R * pos_ptb_vec;                                
                res_ptb = obj.getAS2Res(R,P_ptb,initParams,kappa,L,lane_idx,dir,SubSegIdx);
                p_jac(:,i) = (res_ptb - anchored_res)/eps;
            end
            p_jac = InvMahalanobis(p_jac,cov);

            % 3: Perturbation to Segment initParams
            for i=1:3
                initParams_ptb_vec = zeros(1,3);
                initParams_ptb_vec(i) = eps;
                initParams_ptb = initParams + initParams_ptb_vec;               
                res_ptb = obj.getAS2Res(R,P,initParams_ptb,kappa,L,lane_idx,dir,SubSegIdx);                
                param_jac(:,i) = (res_ptb - anchored_res)/eps;
            end
            param_jac = InvMahalanobis(param_jac,cov);

            % 4: Perturbation to Sub-Segment Params (kappa, L)
            for i=1:SubSegIdx
                kappa_ptb_vec = zeros(1,length(kappa));
                kappa_ptb_vec(i) = eps_kappa;
                kappa_ptb = kappa + kappa_ptb_vec;     
                
                L_ptb_vec = zeros(1,length(L));
                L_ptb_vec(i) = eps_L;
                % L positive constraint
                % actual optimization variable is sqrt(L)
                L_ptb = ((L - minL).^(1/2) + L_ptb_vec).^2 + minL;   

                res_ptbK = obj.getAS2Res(R,P,initParams,kappa_ptb,L,lane_idx,dir,SubSegIdx); 
                res_ptbL = obj.getAS2Res(R,P,initParams,kappa,L_ptb,lane_idx,dir,SubSegIdx); 

                kappa_jac(:,i) = (res_ptbK - anchored_res)/eps_kappa;
                L_jac(:,i) = (res_ptbL - anchored_res)/eps_L;
            end
            kappa_jac = InvMahalanobis(kappa_jac,cov);
            L_jac = InvMahalanobis(L_jac,cov);
            
            anchored_res = InvMahalanobis(anchored_res,cov);
        end
        
        %% Arc Spline based measurement 2 residual Phase 2
        function res = getAS2Res(obj,R,P,initParams,kappa,L,lane_idx,dir,SubSegIdx)
            rpy = dcm2rpy(R); psi = rpy(3);
            R2d = [cos(psi) -sin(psi);
                   sin(psi) cos(psi)];

            % Propgate arc node point until the last sub-segment
%             x = initParams(1); y = initParams(2); heading = initParams(3);   
%             for i=1:SubSegIdx
%                 x = x + 1/kappa(i) * (sin(heading + kappa(i) * L(i)) - sin(heading));
%                 y = y - 1/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));
%                 heading = heading + kappa(i) * L(i);
%             end
            if dir == 1
                dij = obj.lane.ly(lane_idx,1);
%                 zij = obj.lane.lz(lane_idx,1);
            elseif dir == 2
                dij = obj.lane.ry(lane_idx,1);
%                 zij = obj.lane.rz(lane_idx,1);
            end
            
            nodePoints = obj.propNode(initParams,kappa,L);
            X = nodePoints(:,SubSegIdx+1);

%             if strcmp(flag,'init')
%                 X = [initParams(1);initParams(2)];                
%             elseif strcmp(flag,'last')
%                 nodePoints = obj.propNode(initParams,kappa,L);
%                 X = nodePoints(:,end);
%             end
            
            Xb = R2d' * (X - P(1:2)); xb = Xb(1); yb = Xb(2);            
            % xb should be 0 (ideal)
            res = [xb;yb - dij];      
        end
        
        %% Update segment tracker
        function obj = updateTracker(obj)
            obj.opt.seg_tracker = 0;
            num = 0;
            for i=1:length(obj.map.arc_segments)
                num = num + 3 + 2 * length(obj.map.arc_segments{i}.kappa);
                % x0, y0, tau0 of every segment are to be optimized
                % kappa and L of every sub-segment are to be optimized
                obj.opt.seg_tracker = [obj.opt.seg_tracker num];
            end
        end
        
        %% Matrix Singularity Analysis
        function errorAnalysis(obj)
            % Function for finding location of matrix singularity
            % 1. Search for any zero in the diagonals
            % 2. Search for repeated row in information matrix and obtain
            % related variables
            % 3. Search for repeated row in Jacobian measurement matrix and
            % obtain related measurements
            % 
            % Using 2 and 3, find repeated arguments in Jacobian &
            % Information matrix
            % 

            a_diags = diag(obj.opt.info);
            idx = find(a_diags == 0);
            disp(idx)
%             figure(2); spy(obj.opt.info); 
             
            % Check for repeated row in information matrix
            % For faster search, we first trim information matrix for lane
            % parameters only and search for repeated candidates
            % 
            n = length(obj.states);
            params_info = obj.opt.info(16*n+1:end,16*n+1:end);
            
            I = []; J = [];
            for i=1:size(params_info,1)
                for j=i+1:size(params_info,1)
                    d = params_info(i,:) - params_info(j,:);
                    if norm(d) < 1e-8
%                         warning(['param_info row ',num2str(i),' and ',num2str(j),' are almost or totally equal'])
                        I = [I i]; J = [J j];
                    end
                end
            end

            for i=1:length(I)
                posidxI = I(i) + 16*n;
                posidxJ = J(i) + 16*n;
                d = obj.opt.info(posidxI,:) - obj.opt.info(posidxJ,:);
                if norm(d) == 0
                    warning(['Information Matrix row ',num2str(posidxI),' and ',num2str(posidxJ),' are identical'])
                    disp(['Parameter Location: ',num2str(I(i)),', ',num2str(J(i))])
                end
            end
            
            error('Information matrix is singular... stopping optimization')
        end
        
        %% Extract Lane Point Covariance Information
        function obj = extractLaneCov(obj)
            obj.opt.LeftCov = zeros(4*length(obj.states),obj.lane.prev_num);
            obj.opt.RightCov = zeros(4*length(obj.states),obj.lane.prev_num);
            
            disp(['Extracting 2D Lane Point Covariance via Multivariate-Delta Method'])
            for i=1:length(obj.states)
                R_idx = 9*(i-1)+1:9*(i-1)+3;
                P_idx = 9*(i-1)+7:9*(i-1)+9;
                state_cov = [obj.opt.cov(R_idx,R_idx), obj.opt.cov(R_idx,P_idx);
                             obj.opt.cov(P_idx,R_idx), obj.opt.cov(P_idx,P_idx)];
                lane_idx = obj.lane.state_idxs(i);
                
                R = obj.states{i}.R;

                for j=1:obj.lane.prev_num
                    LeftL_cov = diag([0, obj.lane.lystd(lane_idx,j)^2, obj.lane.lzstd(lane_idx,j)^2]);
                    RightL_cov = diag([0, obj.lane.lystd(lane_idx,j)^2, obj.lane.lzstd(lane_idx,j)^2]);
                    LeftL = obj.states{i}.left(:,j);
                    RightL = obj.states{i}.right(:,j);

                    CovL = blkdiag(state_cov,LeftL_cov);
                    CovR = blkdiag(state_cov,RightL_cov);

                    jacL = [-R * skew(LeftL);R;R];
                    jacR = [-R * skew(RightL);R;R];

                    GCovL = jacL' * CovL * jacL;
                    GCovR = jacR' * CovR * jacR;

                    obj.opt.LeftCov(4*i-3:4*i,j) = reshape(GCovL(1:2,1:2),4,1);
                    obj.opt.RightCov(4*i-3:4*i,j) = reshape(GCovR(1:2,1:2),4,1);
                end
            end            
        end

    end
    
    %% Static Methods
    methods (Static)
        %% Coefficient matrices for covariance propagation
        function [Ak, Bk] = getCoeff(delRik, delRkkp1, a_k, w_k, dt_k)
            Ak = zeros(9,9); Bk = zeros(9,6);
            Ak(1:3,1:3) = delRkkp1';
            Ak(4:6,1:3) = -delRik * skew(a_k) * dt_k;
            Ak(4:6,4:6) = eye(3);
            Ak(7:9,1:3) = -1/2 * delRik * skew(a_k) * dt_k^2;
            Ak(7:9,4:6) = eye(3) * dt_k;
            Ak(7:9,7:9) = eye(3);
        
            Bk(1:3,1:3) = RightJac(w_k * dt_k) * dt_k;
            Bk(4:6,4:6) = delRik * dt_k;
            Bk(7:9,4:6) = 1/2 * delRik * dt_k^2;
        end

        %% IMU Jacobians
        function [Ji,Jj,Jb] = getIMUJac(resR,Ri,Rj,Vi,Vj,Pi,Pj,bgdi,JdelRij_bg,JdelVij_bg,JdelVij_ba,JdelPij_bg,JdelPij_ba,dtij)
            grav = [0;0;-9.81];
            Ji = zeros(9,9); Jj = zeros(9,9); Jb = zeros(9,6);
            Ji(1:3,1:3) = -InvRightJac(resR) * Rj' * Ri;
            Ji(4:6,1:3) = skew(Ri' * (Vj - Vi - grav * dtij));
            Ji(4:6,4:6) = -Ri';
            Ji(7:9,1:3) = skew(Ri' * (Pj - Pi - Vi * dtij - 1/2 * grav * dtij^2));
            Ji(7:9,4:6) = -Ri' * dtij;
            Ji(7:9,7:9) = -eye(3);
        
            Jj(1:3,1:3) = InvRightJac(resR);
            Jj(4:6,4:6) = Ri';
            Jj(7:9,7:9) = Ri' * Rj;
        
            Jb(1:3,1:3) = -InvRightJac(resR) * Exp_map(resR)' * RightJac(JdelRij_bg * bgdi) * JdelRij_bg;
            Jb(4:6,1:3) = -JdelVij_bg;
            Jb(4:6,4:6) = -JdelVij_ba;
            Jb(7:9,1:3) = -JdelPij_bg;
            Jb(7:9,4:6) = -JdelPij_ba;
        end
        
        %% WSS Jacobians
        function [Jrv,Jbg,Js] = getWSSJac(Ri,Vi,Si,wi,bgi,bgdi,L)
            Jrv = zeros(3,6); % Horizontally Augmented for R and V
            Jrv(:,1:3) = 1/Si * skew(Ri' * Vi);
            Jrv(:,4:6) = 1/Si * Ri';
        
            Jbg = -1/Si * skew(L);
            Js = -1/Si^2 * (Ri' * Vi - skew(wi - bgi - bgdi) * L);
        end
        
        %% Compute Dog-Leg Step for Approximate Trust Region Method
        function h_dl = ComputeDogLeg(h_gn,h_gd,tr_rad)
            N_hgn = norm(h_gn); N_hgd = norm(h_gd);
            if N_hgn <= tr_rad
                h_dl = h_gn;
            elseif N_hgd >= tr_rad
                h_dl = (tr_rad/N_hgd) * h_gd;
            else
                v = h_gn - h_gd;
                hgdv = h_gd' * v;
                vsq = v' * v;
                beta = (-hgdv + sqrt(hgdv^2 + (tr_rad^2 - h_gd' * h_gd)*vsq)) / vsq;
                h_dl = h_gd + beta * v;
            end
        end
        
        %% Update trust-region radius according to gain ratio value (rho)
        function Delta = UpdateDelta(rho,tr_rad,eta1,eta2,gamma1,gamma2)
            if rho >= eta2
                Delta = gamma2 * tr_rad;
            elseif rho < eta1
                Delta = gamma1 * tr_rad;
            else
                Delta = tr_rad;
            end
        end
        
        %% Update Levenberg-Marquardt lambda value according to gain ration value
        function Lambda = UpdateLambda(lambda,Lu,Ld,flag)
            if flag
                Lambda = max([lambda/Ld,1e-7]);
            else
                Lambda = min([lambda*Lu,1e7]);
            end
        end
        
        %% Detect oscillation by computing maximum deviation from average (recent 5 costs)
        function osc_flag = DetectOsc(cost_stack)
            last_five = cost_stack(end-4:end);
            avg = mean(last_five);
            delta = last_five - avg;
            % If all recent costs are near the average, detect oscillation
            if max(delta) < 1e2
                osc_flag = true;
            else
                osc_flag = false;
            end
        end
        
        %% Find Maximum column index
        function [idx, dir] = maxColIdx(arr,num)
            [r,c] = find(arr==num);
            if isempty(c)
                error('No matched number')
            else
                [idx,max_idx] = max(c);
                if r(max_idx) == 1
                    dir = 'left';
                elseif r(max_idx) == 2
                    dir = 'right';
                end
            end
        end

        %% Find Minimum column index
        function [idx, dir] = minColIdx(arr,num)
            [r,c] = find(arr==num);
            if isempty(c)
                error('No matched number')
            else
                [idx,min_idx] = min(c);
                if r(min_idx) == 1
                    dir = 'left';
                elseif r(min_idx) == 2
                    dir = 'right';
                end
            end
        end
        
        %% Propagate Arc Center Points
        function centerPoints = propCenter(initParams,kappa,L)
            x0 = initParams(1); y0 = initParams(2); heading = initParams(3);
            
            for i=1:length(kappa)
                if i == 1
                    centerPoints = [x0 - 1/kappa(i) * sin(heading); y0 + 1/kappa(i) * cos(heading)];
                else
                    centerPoints = [centerPoints centerPoints(:,end) + (1/kappa(i-1) - 1/kappa(i)) * [sin(heading);-cos(heading)]];
                end
                heading = heading + kappa(i) * L(i);
            end
        end
        
        %% Propagate Arc Node Points
        function nodePoints = propNode(initParams,kappa,L)
            x0 = initParams(1); y0 = initParams(2); heading = initParams(3);            
            nodePoints = [x0; y0];
            
            for i=1:length(kappa)
                nodePoints = [nodePoints nodePoints(:,end) + 1/kappa(i) * [sin(heading + kappa(i) * L(i)) - sin(heading);
                                                                           -cos(heading + kappa(i) * L(i)) + cos(heading)]];
                heading = heading + kappa(i) * L(i);
            end
        end
        
        %% Compute Confidence Ellipse Point
        function pc = ConfidenceEllipse(mu, Sigma, p)
            s = -2 * log(1 - p);
            
            pc = zeros(2,100);            
            
            sig = reshape(Sigma,2,2);
            [V, D] = eig(sig * s);
        
            t = linspace(0, 2 * pi);
            a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        
            pc(1,:) = a(1,:) + mu(1);
            pc(2,:) = a(2,:) + mu(2);             
        end

    end
end

