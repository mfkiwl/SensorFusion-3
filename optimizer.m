classdef optimizer < handle
% OPTIMIZER - Sensor Fusion Optimizer Module
%   All input data are from the dataprocessor module
%   [Usage Format: an example]
%   sol = OPTIMIZER(imu,gnss,lane,can,imu_bias,t,prior_covs,'basic',options)
%
%   [imu, gnss, lane, can, imu_bias, t] : Pre-processed data
% 
%   [Mode] : 'basic', 'partial1', 'partial2', 'full'
%   * Basic: INS + GNSS Fusion
%   * Partial1: INS + GNSS + WSS Fusion
%   * Partial2: INS + GNSS + Lane Fusion
%   * Full: INS + GNSS + WSS + Lane Fusion 
% 
%   [Options] : struct variable containing NLS optimization options
%   * options.CostThres: Cost difference threshold
%   * options.StepThres: Step size threshold
%   * options.IterThres: Iteration limit
%   * options.Algorithm: GN(Gauss-Newton),TR(Trust-Region)
%   ** GN : no guarantee of convergence to local minima
%   When using TR as solver algorithm, need to define parameters.
%   - TR parameters at 'main.m' should work fine
%   - There are cases where TR algorithm converges to a weird local optima
%     --> Tune model covariance values 
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
        
        validL
        lambdaL % Left Lane Lambda Array
        validL_map
        I_L
        J_L
        validR
        lambdaR % Right Lane Lambda Array
        validR_map
        I_R
        J_R

        bias = {} % bias variables saver
        bgd_thres = 1e-2; % Gyro bias repropagation threshold
        bad_thres = 1e-1; % Accel bias repropagation threshold
        
%         map % map variable for 'full' or '2-phase' mode
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
            
            n = length(obj.imu)+1;
            for i=1:n
                bias = struct();
                bias.bg = obj.imu_bias.gyro';
                bias.ba = obj.imu_bias.accel';
                bias.bgd = zeros(3,1);
                bias.bad = zeros(3,1);

                obj.bias = [obj.bias {bias}];
            end
            
            % Vehicle Lever Arm Vector
            obj.params =  [1.5; 0; 0.6]; 
            
            % Lambda Array for lane matching
            obj.lambdaL = zeros(n,obj.lane.prev_num);
            obj.validL_map = zeros(n,obj.lane.prev_num);
            obj.lambdaR = zeros(n,obj.lane.prev_num);
            obj.validR_map = zeros(n,obj.lane.prev_num);
            
            % Validate Lane Measurement based on std value and probability
            obj.validateLane();

            % Perform Pre-integration
            obj.integrate(1:n-1);

            % Perform INS propagation
            obj.ins();
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
%                 obj.map.optimize();
% 
%                 obj.map.segment_info = zeros(2,size(obj.lane.FactorValidIntvs,1));
% 
%                 LeftSegNum = 1; RightSegNum = 2;
%                 obj.map.segment_info(1,1) = LeftSegNum;
%                 obj.map.segment_info(2,1) = RightSegNum;
%     
%                 for i=1:length(obj.lane.LC_dirs)
%                     if strcmp(obj.lane.LC_dirs{i},'left')
%                         tmp = LeftSegNum;
%                         LeftSegNum = max([LeftSegNum,RightSegNum]) + 1;
%                         RightSegNum = tmp;
%                     elseif strcmp(obj.lane.LC_dirs{i},'right')
%                         tmp = RightSegNum;
%                         RightSegNum = max([LeftSegNum, RightSegNum]) + 1;
%                         LeftSegNum = tmp;
%                     end
%                     obj.map.segment_info(1,i+1) = LeftSegNum;
%                     obj.map.segment_info(2,i+1) = RightSegNum;            
%                 end
            end           
        end

        %% Optimize
        function obj = optimize(obj)
            % Solve Sparse Nonlinear Least Squares in Batch Manner                     

            disp('[Optimization Starts...]')

            n = length(obj.states);

            
            % Run optimization depending on algorithm options
            % [Basic Mode]: INS + GNSS Fusion
            % [Partial Mode]: INS + GNSS + WSS Fusion
            % [Partial Mode2]: INS + GNSS + Lane Fusion
            % [Full Mode]: INS + GNSS + WSS + Lane Fusion
            
            if strcmp(obj.mode,'basic')
                obj.opt.x0 = zeros(15*n,1);
            elseif strcmp(obj.mode,'partial1')
                obj.opt.x0 = zeros(16*n+3,1);   
            elseif strcmp(obj.mode,'partial2')
                obj.opt.x0 = zeros(15*n + length(obj.I_L) + length(obj.I_R),1);  
            elseif strcmp(obj.mode,'full')
                % Test if lever arm vector can be estimated
                obj.opt.x0 = zeros(16*n + length(obj.I_L) + length(obj.I_R) + 3,1); 
            end

            if strcmp(obj.opt.options.Algorithm,'GN')
                obj.GaussNewton();
            elseif strcmp(obj.opt.options.Algorithm,'TR')
                obj.TrustRegion();
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

%                 if strcmp(obj.mode,'full')
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
                plot(left_(1,1),left_(2,1),'c.');
                plot(right_(1,1),right_(2,1),'m.');
            end

            for i=1:size(obj.lane.FactorValidIntvs,1)
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                plot(obj.states{lb}.P(1),obj.states{lb}.P(2),'gs')
                plot(obj.states{ub}.P(1),obj.states{ub}.P(2),'gs');
            end

%             if strcmp(obj.mode,'full')
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
            p_est = geoplot(P_lla(:,1),P_lla(:,2),'g.'); 
            
%             if strcmp(obj.mode,'full')
%                 n = length(obj.map.arc_segments);
%                 for i=1:n
%                     m = length(obj.map.arc_segments{i}.kappa);
%                     seg = obj.map.arc_segments{i};
%                     heading = seg.tau0;
%                     SegPoints = [seg.x0;
%                                  seg.y0];
%                     for j=1:m
%                         kappa = seg.kappa(j); L = seg.L(j);
%                         headingPrev = heading;
%                         heading = heading + kappa * L;
%                         headingCurr = heading;
%     
%                         heading_ = linspace(headingPrev,headingCurr,1e3);
%                         addedSegPoints = SegPoints(:,end) + 1/kappa * [sin(heading_) - sin(headingPrev);
%                                                                        -cos(heading_) + cos(headingPrev)];
%                         SegPoints = [SegPoints addedSegPoints];
%                         
%                         lla_node1 = enu2lla([addedSegPoints(1,1),addedSegPoints(2,1),0],obj.gnss.lla0,'ellipsoid');
%                         lla_node2 = enu2lla([addedSegPoints(1,end),addedSegPoints(2,end),0],obj.gnss.lla0,'ellipsoid');
%                         p_node = geoplot(lla_node1(1),lla_node1(2),'co');
%                         geoplot(lla_node2(1),lla_node2(2),'co');
%                     end
%                     
%                     lla_seg = enu2lla([SegPoints' zeros(size(SegPoints,2),1)],obj.gnss.lla0,'ellipsoid');
%                     lla_node3 = enu2lla([SegPoints(1,1),SegPoints(2,1),0],obj.gnss.lla0,'ellipsoid');
%                     lla_node4 = enu2lla([SegPoints(1,end),SegPoints(2,end),0],obj.gnss.lla0,'ellipsoid');
%     
%                     p_lane = geoplot(lla_seg(:,1),lla_seg(:,2),'k-');
%                     p_seg = geoplot(lla_node3(1),lla_node3(2),'ms');
%                     geoplot(lla_node4(1),lla_node4(2),'ms');
%                     
%                     
%                 end
%                 geobasemap satellite
% 
%                 title('Optimized Trajectory Comparison')
%                 legend([p_est,p_gnss,p_lane,p_seg,p_node],...
%                        'Optimized Vehicle Trajectory','GNSS Measurements',...
%                        'Continuous Piecewise Arc Spline',...
%                        'Arc Segment Node','Arc Sub-Segment Node')
% %                 title('Scenario Trajectory')
% %                 legend([p_gnss],'GNSS Measurements')
%             else
%                 geobasemap satellite
% 
%                 title('Optimized Trajectory Comparison')
%                 legend([p_est,p_gnss],...
%                        'Optimized Vehicle Trajectory','GNSS Measurements')
%             end
            geobasemap satellite

            title('Optimized Trajectory Comparison')
            legend([p_est,p_gnss],...
                   'Optimized Vehicle Trajectory','GNSS Measurements')

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
        function visualizeHD(obj,varargin)
            figure(20); hold on; grid on; axis equal;
            Lane1 = varargin{1}; % shp file 1
            Lane2 = varargin{2}; % shp file 2
            
%             roadspec = makesymbolspec('Line',...
%                                        {'Default','Color','blue','LineStyle','--','LineWidth',1}, ...
%                                        {'Kind',"503",'Color','green','LineStyle','--','LineWidth',1});

            p_l1 = mapshow(Lane1,'Color','blue','LineStyle','--','Marker','.','MarkerEdgeColor','red'); 
%             p_l1 = mapshow(Lane1,'SymbolSpec',roadspec);
            p_l2 = mapshow(Lane2,'Color','blue','LineStyle','--','Marker','.','MarkerEdgeColor','red');
            % Convert vehicle states to UTM-K coords            
            
            ENU = [];
            for i=1:length(obj.states)
                ENU = [ENU; obj.states{i}.P'];                
            end
            LLA = enu2lla(ENU,obj.gnss.lla0,'ellipsoid');
            p = projcrs(32652); % Projected Coordinate Reference System
            [x,y] = projfwd(p,LLA(:,1),LLA(:,2));
            p_v = plot(x,y,'k.');
            
            
            if length(varargin) > 2 % ArcMap is given
                for i=1:length(varargin)-2
                    ParamMap = varargin{i+2};
                    for j=1:length(ParamMap.arc_segments)
                        seg = ParamMap.arc_segments{j}.params;
                        enu_init = [seg.x0, seg.y0, 0];
                        lla_init = enu2lla(enu_init,obj.gnss.lla0,'ellipsoid');
                        [x_init,y_init] = projfwd(p,lla_init(1),lla_init(2)); % Initial Point Transferred
                        heading = seg.tau0;
                        SegPoints = [seg.x0; seg.y0];
                        SegPointsUTMK = [x_init;y_init];
                        for k=1:length(seg.kappa)
                            kappa = seg.kappa(k); L = seg.L(k);
                            headingPrev = heading;
                            heading = heading + kappa * L;
                            headingCurr = heading;

                            heading_ = linspace(headingPrev,headingCurr,1e3);
                            addedSegPoints = SegPoints(:,end) + 1/kappa * [sin(heading_) - sin(headingPrev);
                                                                           -cos(heading_) + cos(headingPrev)];
                            
                            % Convert addedSegPoints : enu --> lla --> utm-k
                            addedSegPointsLLA = enu2lla([addedSegPoints;zeros(1,size(addedSegPoints,2))]',obj.gnss.lla0,'ellipsoid');
                            [addedSegPointsUTMKX,addedSegPointsUTMKY] = projfwd(p,addedSegPointsLLA(:,1),addedSegPointsLLA(:,2));

                            SegPoints = [SegPoints, addedSegPoints];
                            SegPointsUTMK = [SegPointsUTMK, [addedSegPointsUTMKX';addedSegPointsUTMKY']];
                            p_node = plot(addedSegPointsUTMKX(1),addedSegPointsUTMKY(1),'co');
                            plot(addedSegPointsUTMKX(end),addedSegPointsUTMKY(end),'co');
                        end
                        p_lane = plot(SegPointsUTMK(1,:),SegPointsUTMK(2,:),'Color',[0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',2);
                        p_seg = plot(SegPointsUTMK(1,1),SegPointsUTMK(2,1),'ms');
                        plot(SegPointsUTMK(1,end),SegPointsUTMK(2,end),'ms');

                    end
                end
                xlabel('UTM-K X (m)'); ylabel('UTM-K Y (m)'); title('Comparison with HD Map')
                legend([p_v,p_l1,p_lane,p_node], ...
                        'Vehicle Trajectory', ...
                        'HD Map Link(Blue)/Node(Red)', ...
                        'Arc Spline Ego Lane', ...
                        'Arc Spline Node Points');
            else
                xlabel('UTM-K X (m)'); ylabel('UTM-K Y (m)'); title('HD Map')
%                 legend([p_l, p_v],'HD Map Lanes','Vehicle Trajectory')
                legend(p_l1,'HD Map Link(Blue)/Node(Red)')
            end


            

            
        end        
        
        %% Visualize with Reference and Comparison * Only when there is reference data
        function visualizeRef(obj,varargin)

            ref = varargin{1};
            if length(varargin) == 5
                vins_stereo_imu = varargin{2};
                vins_stereo = varargin{3};
                vins_mono_imu = varargin{4};
                dso = varargin{5};
            elseif length(varargin) == 4
                vins_stereo_imu = varargin{2};
                vins_stereo = varargin{3};
                vins_mono_imu = varargin{4};
                dso = [];
            else
                error('Invalid input size')
            end

            % Number of position data should be equal to number of states
            obj.opt.ENU_ref = lla2enu(ref.pos,obj.gnss.lla0,'ellipsoid');
            n = length(obj.states);
            if size(obj.opt.ENU_ref,1) ~= n
                error('Data size should be matched!')
            end

            obj.opt.E = zeros(1,n);
            for i=1:n
                P = obj.states{i}.P;
                obj.opt.E(i) = norm(obj.opt.ENU_ref(i,1:2)' - P(1:2));
            end
            obj.opt.SE =obj.opt.E.^2;
            obj.opt.MSE = mean(obj.opt.SE);
            obj.opt.RMSE = sqrt(obj.opt.MSE);
            
            disp('=========Optimization Summary with Reference data=========')
            if strcmp(obj.mode,'basic')
                disp('vSensor Fusion Mode (INS + GNSS)----')
            elseif strcmp(obj.mode,'partial')
                disp('----Sensor Fusion Mode (INS + GNSS + WSS)----')
            end
            
            disp(['RMSE: ',num2str(obj.opt.RMSE,3),'m'])
            disp(['Max Error: ',num2str(max(obj.opt.E),3),'m'])
            

            P = [];

            for i=1:n
                P_ = obj.states{i}.P;
                P = [P P_];                
            end
            figure(30);
            p_gnss = geoplot(obj.gnss.pos(:,1),obj.gnss.pos(:,2),'b.'); grid on; hold on;
            P_lla = enu2lla(P',obj.gnss.lla0,'ellipsoid');                        
            obj.opt.P_lla = P_lla;
            p_est = geoplot(P_lla(:,1),P_lla(:,2),'r.'); 
            p_ref = geoplot(ref.pos(:,1),ref.pos(:,2),'k.');
            geobasemap satellite

            title('Optimized Trajectory Comparison')
            legend([p_est,p_gnss,p_ref],...
                   'Optimized Vehicle Trajectory','GNSS Measurements','Ground Truth')
            
            
            % Plot on ENU frame
            figure(31); grid on; hold on; axis equal;
            p_ref = plot(obj.opt.ENU_ref(:,1),obj.opt.ENU_ref(:,2),'k-.','Marker','o','MarkerSize',2);
            p_est = plot(P(1,:),P(2,:),'r--');
            p_vins_st_imu = plot(vins_stereo_imu.x,vins_stereo_imu.y,'m--');
            p_vins_st = plot(vins_stereo.x,vins_stereo.y,'c--');
            p_vins_mono = plot(vins_mono_imu.x,vins_mono_imu.y,'g--');
            
            xlabel('Global X(m)','FontSize',14); ylabel('Global Y(m)','FontSize',14);
            title('Vehicle Trajectory Reconstruction','FontSize',14);
            if ~isempty(dso)
                p_dso = plot(dso.x,dso.y,'b--');
                lgd = legend([p_ref,p_est,p_vins_st_imu,p_vins_st,p_vins_mono,p_dso], ...
                             'Reference Trajectory', ...
                             'Proposed Sensor Fusion', ...
                             'VINS Stereo + IMU', ...
                             'VINS Stereo','VINS Mono + IMU','DSO');
                lgd.FontSize = 14;
            else
                lgd = legend([p_ref,p_est,p_vins_st_imu,p_vins_st,p_vins_mono], ...
                             'Reference Trajectory', ...
                             'Proposed Sensor Fusion', ...
                             'VINS Stereo + IMU', ...
                             'VINS Stereo','VINS Mono + IMU');
                lgd.FontSize = 14;
            end
            
            
            % Error Analysis for other algorithms
            obj.opt.comparison_error = struct();
            % VINS Stereo IMU    
            obj.opt.comparison_error.vins_stimu = obj.computeError(obj.opt.ENU_ref(:,1:2)',[vins_stereo_imu.x'; vins_stereo_imu.y']);
            SE = obj.opt.comparison_error.vins_stimu.^2;
            MSE = mean(SE); RMSE = sqrt(MSE);
            disp('----VINS Stereo + IMU----')
            disp(['RMSE: ',num2str(RMSE,3),'m'])
            disp(['Max Error: ',num2str(max(obj.opt.comparison_error.vins_stimu),3),'m'])

            % VINS Stereo 
            obj.opt.comparison_error.vins_st = obj.computeError(obj.opt.ENU_ref(:,1:2)',[vins_stereo.x'; vins_stereo.y']);
            SE = obj.opt.comparison_error.vins_st.^2;
            MSE = mean(SE); RMSE = sqrt(MSE);
            disp('----VINS Stereo----')
            disp(['RMSE: ',num2str(RMSE,3),'m'])
            disp(['Max Error: ',num2str(max(obj.opt.comparison_error.vins_st),3),'m'])
            
            % VINS Mono
            obj.opt.comparison_error.vins_mnimu = obj.computeError(obj.opt.ENU_ref(:,1:2)',[vins_mono_imu.x'; vins_mono_imu.y']);
            SE = obj.opt.comparison_error.vins_mnimu.^2;
            MSE = mean(SE); RMSE = sqrt(MSE);
            disp('----VINS Mono + IMU----')
            disp(['RMSE: ',num2str(RMSE,3),'m'])
            disp(['Max Error: ',num2str(max(obj.opt.comparison_error.vins_mnimu),3),'m'])

            % DSO (If existing)
            if ~isempty(dso)
                obj.opt.comparison_error.dso = obj.computeError(obj.opt.ENU_ref(:,1:2)',[dso.x'; dso.y']);
                SE = obj.opt.comparison_error.dso.^2;
                MSE = mean(SE); RMSE = sqrt(MSE);
                disp('----DSO----')
                disp(['RMSE: ',num2str(RMSE,3),'m'])
                disp(['Max Error: ',num2str(max(obj.opt.comparison_error.dso),3),'m'])
            else
                disp('----DSO----')
                disp('Algorithm Failed')
            end

        end

        %% Plot lane point std values
        function plotConfEllipse(obj,p)
            % 0 < p < 1 : Confidence Level
            % Left Lane
            figure(200); grid on; axis equal; hold on;
            % i=1:length(obj.states) --> too much data! 
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
        %% Label Valid Lane Indices
        function obj = validateLane(obj)
            % Initialize lambda only for valid lane measurements
            % All labeled lane variables are to be optimized for 'partial2'
            % or 'full' modes.
            n = length(obj.imu)+1;
            for i=1:n-1
%                 disp(i)
                laneIdx1 = obj.lane.state_idxs(i);
                laneIdx2 = obj.lane.state_idxs(i+1);
                
                % Left Lane
                % Skipped if lane reliability is too low
                if obj.lane.prob(laneIdx1,2) > obj.lane.prob_thres && obj.lane.prob(laneIdx2,2) > obj.lane.prob_thres

                    for j=1:obj.lane.prev_num-1
                        Qij_std = obj.lane.lystd(laneIdx1,j);
                        Qijp1_std = obj.lane.lystd(laneIdx1,j+1);
                        Qip1j_std = obj.lane.lystd(laneIdx2,j);

                        % Valid if measurements for Qij, Qijp1, Qip1j have
                        % low enough std value
                        if Qij_std < obj.lane.std_thres && Qijp1_std < obj.lane.std_thres && Qip1j_std < obj.lane.std_thres
                            obj.validL(i,j) = 1;
                            obj.validL(i,j+1) = 1;
                            obj.validL(i+1,j) = 1;
                            obj.lambdaL(i,j) = 1;
                        end                        
                    end                
                end
                
                % Right Lane
                % Skipped if lane reliability is too low
                if obj.lane.prob(laneIdx1,3) > obj.lane.prob_thres && obj.lane.prob(laneIdx2,3) > obj.lane.prob_thres
                    
                    for j=1:obj.lane.prev_num-1
                        Qij_std = obj.lane.rystd(laneIdx1,j);
                        Qijp1_std = obj.lane.rystd(laneIdx1,j+1);
                        Qip1j_std = obj.lane.rystd(laneIdx2,j);
                        
                        % Valid if measurements for Qij, Qijp1, Qip1j have
                        % low enough std value
                        if Qij_std < obj.lane.std_thres && Qijp1_std < obj.lane.std_thres && Qip1j_std < obj.lane.std_thres
                            obj.validR(i,j) = 1;
                            obj.validR(i,j+1) = 1;
                            obj.validR(i+1,j) = 1;
                            obj.lambdaR(i,j) = 1;
                        end                        
                    end                
                end
            end

            [obj.I_L,obj.J_L] = find(obj.validL);
            [obj.I_R,obj.J_R] = find(obj.validR);

            % Map linear indices
            % Left Lane
            cntL = 1;
            for i=1:length(obj.I_L)
                I = obj.I_L(i); J = obj.J_L(i);
                obj.validL_map(I,J) = cntL;
                cntL = cntL + 1;
            end

            % Right Lane
            cntR = 1;
            for i=1:length(obj.I_R)
                I = obj.I_R(i); J = obj.J_R(i);
                obj.validR_map(I,J) = cntR;
                cntR = cntR + 1;
            end
        end

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
        %
        % * Vehicle Body Frame: [x, y, z] = [Forward, Left, Up]
        % * World Frame: [x, y, z] = [East, North, Up]
        % * Camera Frame(Phone): https://source.android.com/docs/core/sensors/sensor-types
        % Note that IMU (Android) measurements have gravitational effects,
        % which should be considered when modeling vehicle 3D kinematics 

        % # Vehicle State Variables
        % # R: Body-to-world frame rotational matrix
        % # V: World frame velocity vector
        % # P: World frame position vector (local frame coords, not geodetic)
        %
        % Implemented by JinHwan Jeon, 2022

            disp('[INS Propagation...]')
            th = pi/2 - pi/180 * obj.gnss.bearing(1);
            R = [cos(th) -sin(th) 0;
                 sin(th) cos(th)  0;
                 0       0        1];
            Vned = obj.gnss.vNED(1,:);
            V = [Vned(2); Vned(1); -Vned(3)];
            if isempty(obj.snap)
                P = lla2enu(obj.gnss.pos(1,:),obj.gnss.lla0,'ellipsoid')';
            else
                P = lla2enu(obj.snap,obj.gnss.lla0,'ellipsoid');
                P = P';
            end
            
            grav = [0;0;-9.81];

            state_ = struct();
            state_.R = R; state_.V = V; state_.P = P;
            state_.WSF = 1;
            
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
            state_.L = obj.params;            
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
                state_.WSF = 1;    
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
                
                state_.lambdaL = zeros(1,prev_num);
                state_.lambdaR = zeros(1,prev_num);

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
%             obj.opt.info = jac' * jac;
%             disp('[Computing Sparse Inverse of Information Matrix...]')
%             obj.opt.cov = sparseinv(obj.opt.info);
%             obj.extractLaneCov();
        end

        %% Trust Region Method
        function cost = TrustRegion(obj)
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
%                 dummy_lambdaL = obj.lambdaL;
%                 dummy_lambdaR = obj.lambdaR;
                dummy_params = obj.params;
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
                    obj.params = dummy_params;
%                     obj.lambdaL = dummy_lambdaL;
%                     obj.lambdaR = dummy_lambdaR;

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
                
                if tr_rad / step_size > 1e10
                    % If Trust Region size is too large compared to
                    % iterative step size, manually reduce trust region
                    % radius to the current step size for efficient
                    % rejection in the future iterations
                    tr_rad = step_size;
                end
            end

%             obj.opt.info = jac' * jac;
%             disp('[Computing Sparse Inverse of Information Matrix...]')
%             obj.opt.cov = sparseinv(obj.opt.info);
%             obj.extractLaneCov();
        end

        %% Optimization Cost Function
        function [res,jac] = cost_func(obj,x0)
            obj.retract(x0,'normal');
            [Pr_res,Pr_jac] = obj.CreatePrBlock();
            [MM_res,MM_jac] = obj.CreateMMBlock();
            [GNSS_res,GNSS_jac] = obj.CreateGNSSBlock();
            [WSS_res,WSS_jac] = obj.CreateWSSBlock();
            [LM_res,LM_jac] = obj.CreateLMBlock();
            
            res = vertcat(Pr_res,MM_res,GNSS_res,WSS_res,LM_res);
            jac = vertcat(Pr_jac,MM_jac,GNSS_jac,WSS_jac,LM_jac);   
            obj.opt.jac = jac;
        end
        
        %% Prior Residual and Jacobian
        function [Pr_res,Pr_jac] = CreatePrBlock(obj)
            % Prior Model
            % 
            % Prior for anchoring initial states and wsf values
            % Implemented by JinHwan Jeon, 2022

            n = length(obj.states);
            blk_width = length(obj.opt.x0);
            blk_height = 15;
            adder = 0;

            if strcmp(obj.mode,'partial1') || strcmp(obj.mode,'full')
                adder = n ;                
                blk_height = blk_height + n;                
            end

            Veh_cov = blkdiag(obj.covs.prior.R,obj.covs.prior.V,obj.covs.prior.P);
            Bias_cov = blkdiag(obj.covs.prior.Bg,obj.covs.prior.Ba);
            
            I = zeros(1,9*3 + 6 + adder); J = I; V = I;

            % Vehicle States Prior
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
            
            % Wheel Scaling Factor Prior
            if ~strcmp(obj.mode,'partial2')
                WSS_res = zeros(n,1);
                for i=1:n
                    WSF = obj.states{i}.WSF;
                    WSS_res(i) = InvMahalanobis(WSF - 1,obj.covs.prior.WSF);
                    I(9*3+6+i) = 15+i;
                    J(9*3+6+i) = 15*n+i;
                    V(9*3+6+i) = InvMahalanobis(1,obj.covs.prior.WSF);
                end
                
                % Lever arm prior
                lin_point = [1.5;0;0.7];
                cov = obj.covs.prior.Params;
                Lev_res = InvMahalanobis(obj.params - lin_point,cov);
                I_Lev = [1,2,3]; J_Lev = [1,2,3]; V_Lev = diag(InvMahalanobis(eye(3),cov))';

                if strcmp(obj.mode,'partial1')
                    Lev_jac = horzcat(sparse(3,16*n),sparse(I_Lev,J_Lev,V_Lev,3,3));
                elseif strcmp(obj.mode,'full')
                    Lev_jac = horzcat(sparse(3,16*n+length(obj.I_L)+length(obj.I_R)),sparse(I_Lev,J_Lev,V_Lev,3,3));
                end                
            else
                WSS_res = []; Lev_res = []; Lev_jac = [];
            end
            
            % Lane Measurement Prior
            if ~strcmp(obj.mode,'partial1')
                if strcmp(obj.mode,'partial2')
                    mul = 15; adder2 = 0;
                elseif strcmp(obj.mode,'full')
                    mul = 16; adder2 = 3;
                end
                % Left Lane
                n_L = length(obj.I_L);
                Lane_resL = zeros(n_L,1);
                % normalize covariance
                I_L_ = zeros(1,n_L); J_L_ = I_L_; V_L_ = I_L_;
                for i=1:n_L
                    I_ = obj.I_L(i); J_ = obj.J_L(i);
                    laneIdx = obj.lane.state_idxs(I_);
                    cov = obj.lane.lystd(laneIdx,J_)^2; 
                    
                    Lane_resL(i) = InvMahalanobis(obj.states{I_}.left(2,J_) - obj.lane.ly(laneIdx,J_),cov);
                    I_L_(i) = i; J_L_(i) = i; 
                    V_L_(i) = InvMahalanobis(1,cov);
                end
                % Right Lane
                n_R = length(obj.I_R);
                Lane_resR = zeros(n_R,1);
                I_R_ = zeros(1,n_R); J_R_ = I_R_; V_R_ = I_R_;
                for i=1:n_R
                    I_ = obj.I_R(i); J_ = obj.J_R(i);
                    laneIdx = obj.lane.state_idxs(I_);
                    cov = obj.lane.rystd(laneIdx,J_)^2;
                    Lane_resR(i) = InvMahalanobis(obj.states{I_}.right(2,J_) - obj.lane.ry(laneIdx,J_),cov);
                    I_R_(i) = i; J_R_(i) = i; 
                    V_R_(i) = InvMahalanobis(1,cov);
                end
                Lane_jacL = horzcat(sparse(n_L,mul*n),sparse(I_L_,J_L_,V_L_,n_L,n_L+n_R+adder2));
                Lane_jacR = horzcat(sparse(n_R,mul*n+n_L),sparse(I_R_,J_R_,V_R_,n_R,n_R+adder2));
                Lane_res = vertcat(Lane_resL,Lane_resR);
                Lane_jac = vertcat(Lane_jacL,Lane_jacR);
            else
                Lane_res = []; Lane_jac = [];
            end            
            Pr_res = vertcat(Veh_res,Bias_res,WSS_res,Lane_res,Lev_res);
            Pr_jac = vertcat(sparse(I,J,V,blk_height,blk_width),Lane_jac,Lev_jac);
        end
        
        %% IMU Residual and Jacobian
        function [MM_res,MM_jac] = CreateMMBlock(obj)
            % IMU Measurment Model
            %
            % Compare relative motion with preintegrated IMU readings 
            % Implemented "On-Manifold Preintegration for Real-Time Visual-Inertial Odometry" 
            % by C. Forster etal,
            % 
            % Implemented in MATLAB by JinHwan Jeon, 2022

            n = length(obj.states);
            blk_width = length(obj.opt.x0);
%             blk_height = 15 * (n-1);
            grav = [0;0;-9.81];

            if strcmp(obj.mode,'partial1') || strcmp(obj.mode,'full')
                I_s = zeros(1,2*(n-1)); J_s = I_s; V_s = I_s;
                s_cov = obj.covs.imu.ScaleFactorNoise;
                WSF_res = zeros(n-1,1);
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
                
                if strcmp(obj.mode,'partial1') || strcmp(obj.mode,'full')
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
            
            if strcmp(obj.mode,'partial1') || strcmp(obj.mode,'full')
                % Augment WSF residual and jacobian
                WSF_jac = sparse(I_s,J_s,V_s,n-1,blk_width);
                
                MM_res = [MM_res; WSF_res];
                MM_jac = [MM_jac; WSF_jac];
            end

        end

        %% GNSS Residual and Jacobian
        function [GNSS_res,GNSS_jac] = CreateGNSSBlock(obj)
            % GNSS Measurement Model
            % 
            % Supports Loosely Coupled GNSS/INS
            % Convert LLA format readings to local ENU and compare with
            % vehicle states
            % 
            % Implemented by JinHwan Jeon, 2022

            idxs = obj.gnss.state_idxs;
            n = length(idxs);
            blk_width = length(obj.opt.x0);
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
            % Wheel Speed Sensor Measurement Model
            %
            % This model is based on the Vehicle Non-holonomic Constraint
            % "Vertical and lateral velocity vector at the rear axis is 0".
            % Based on the 3D Rigid Body Kinematic equation, residual
            % function was designed.
            % 
            % Implemented by JinHwan Jeon, 2022

            if strcmp(obj.mode,'basic') || strcmp(obj.mode,'partial2')
                WSS_res = [];
                WSS_jac = [];
            else
                n = length(obj.states);
                blk_width = length(obj.opt.x0);
                if strcmp(obj.mode,'partial1')
                    var_idx = 16*n;
                elseif strcmp(obj.mode,'full')
                    var_idx = 16*n+length(obj.I_L)+length(obj.I_R);
                end
                L = obj.params;
    
                WSS_res = zeros(3*(n-2),1);
%                 m = 9*3+3;
                m = 9*4+3;
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
    
                    [Jrv,Jbg,Js,JL] = obj.getWSSJac(Ri,Vi,Si,wi,bgi,bgdi,L);

                    [I_rv,J_rv,V_rv] = sparseFormat(3*(i-2)+1:3*(i-2)+3,9*(i-1)+1:9*(i-1)+6,InvMahalanobis(Jrv,wss_cov));
                    [I_bg,J_bg,V_bg] = sparseFormat(3*(i-2)+1:3*(i-2)+3,9*n+6*(i-1)+1:9*n+6*(i-1)+3,InvMahalanobis(Jbg,wss_cov));
                    [I_s,J_s,V_s] = sparseFormat(3*(i-2)+1:3*(i-2)+3,15*n+i,InvMahalanobis(Js,wss_cov));
                    [I_Lev,J_Lev,V_Lev] = sparseFormat(3*(i-2)+1:3*(i-2)+3,var_idx+1:var_idx+3,InvMahalanobis(JL,wss_cov));                    
                    
                    I(m*(i-2)+1:m*(i-1)) = [I_rv, I_bg, I_s, I_Lev];
                    J(m*(i-2)+1:m*(i-1)) = [J_rv, J_bg, J_s, J_Lev];
                    V(m*(i-2)+1:m*(i-1)) = [V_rv, V_bg, V_s, V_Lev];
                end

                WSS_jac = sparse(I,J,V,3*(n-2),blk_width);
            end
        end
        
        %% Lane Residual and Jacobian
        function [LM_res,LM_jac] = CreateLMBlock(obj)
            % 3D Lane Factor 
            if strcmp(obj.mode,'basic') || strcmp(obj.mode,'partial1')
                LM_res = []; LM_jac = [];
            else
                if strcmp(obj.mode,'partial2')
                    mul = 15;
                elseif strcmp(obj.mode,'full')
                    mul = 16;
                end
                blk_heightL = nnz(obj.lambdaL);
                blk_heightR = nnz(obj.lambdaR);
                blk_width = length(obj.opt.x0);
                
                [I_L_,J_L_] = find(obj.lambdaL);
                [I_R_,J_R_] = find(obj.lambdaR);

                n_L = length(I_L_); n_R = length(I_R_);
                LM_resL = zeros(blk_heightL,1); LM_resR = zeros(blk_heightR,1);
                I_JacL = zeros(15*n_L,1); J_JacL = I_JacL; V_JacL = I_JacL;
                I_JacR = zeros(15*n_R,1); J_JacR = I_JacR; V_JacR = I_JacR;
                % Left Lane
                for i=1:n_L
                    I = I_L_(i); J = J_L_(i);                    
                    lane_idx = obj.lane.state_idxs(I+1);

                    Ri = obj.states{I}.R; Pi = obj.states{I}.P;
                    Lij = obj.states{I}.left(:,J); 
                    var_idx1 = obj.validL_map(I,J);

                    Lijp1 = obj.states{I}.left(:,J+1);
                    var_idx2 = obj.validL_map(I,J+1);

                    Rip1 = obj.states{I+1}.R; Pip1 = obj.states{I+1}.P;
                    Lip1j = obj.states{I+1}.left(:,J);
                    var_idx3 = obj.validL_map(I+1,J);

                    % Get Residual and Jacobian
                    [res,Jri,Jpi,Jrip1,Jpip1,Jlij,Jlijp1,Jlip1j] = obj.getLaneJac(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J);
                    cov = obj.lane.lystd(lane_idx,J)^2;
                    
                    % Residual
                    LM_resL(i) = InvMahalanobis(res,cov);
                    
                    % Jacobian
                    [I_L1,J_L1,V_L1] = sparseFormat(i,9*(I-1)+1:9*(I-1)+3,InvMahalanobis(Jri,cov));
                    [I_L2,J_L2,V_L2] = sparseFormat(i,9*(I-1)+7:9*(I-1)+9,InvMahalanobis(Jpi,cov));
                    [I_L3,J_L3,V_L3] = sparseFormat(i,9*(I)+1:9*(I)+3,InvMahalanobis(Jrip1,cov));
                    [I_L4,J_L4,V_L4] = sparseFormat(i,9*(I)+7:9*(I)+9,InvMahalanobis(Jpip1,cov));
                    [I_L5,J_L5,V_L5] = sparseFormat(i,mul*length(obj.states)+var_idx1,InvMahalanobis(Jlij,cov));
                    [I_L6,J_L6,V_L6] = sparseFormat(i,mul*length(obj.states)+var_idx2,InvMahalanobis(Jlijp1,cov));
                    [I_L7,J_L7,V_L7] = sparseFormat(i,mul*length(obj.states)+var_idx3,InvMahalanobis(Jlip1j,cov));
                    I_JacL(15*(i-1)+1:15*i) = [I_L1,I_L2,I_L3,I_L4,I_L5,I_L6,I_L7];
                    J_JacL(15*(i-1)+1:15*i) = [J_L1,J_L2,J_L3,J_L4,J_L5,J_L6,J_L7];
                    V_JacL(15*(i-1)+1:15*i) = [V_L1,V_L2,V_L3,V_L4,V_L5,V_L6,V_L7];
                end
                
                LM_jacL = sparse(I_JacL,J_JacL,V_JacL,blk_heightL,blk_width);

                % Right Lane
                for i=1:n_R
                    I = I_R_(i); J = J_R_(i);
                    lane_idx = obj.lane.state_idxs(I+1);

                    Ri = obj.states{I}.R; Pi = obj.states{I}.P;
                    Lij = obj.states{I}.right(:,J);
                    var_idx1 = obj.validR_map(I,J)+length(obj.I_L);

                    Lijp1 = obj.states{I}.right(:,J+1);
                    var_idx2 = obj.validR_map(I,J+1)+length(obj.I_L);

                    Rip1 = obj.states{I+1}.R; Pip1 = obj.states{I+1}.P;
                    Lip1j = obj.states{I+1}.right(:,J);
                    var_idx3 = obj.validR_map(I+1,J)+length(obj.I_L);

                    % Get Residual and Jacobian
                    [res,Jri,Jpi,Jrip1,Jpip1,Jlij,Jlijp1,Jlip1j] = obj.getLaneJac(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J);
                    cov = obj.lane.rystd(lane_idx,J)^2;
                    
                    % Residual
                    LM_resR(i) = InvMahalanobis(res,cov);
                    
                    % Jacobian
                    [I_R1,J_R1,V_R1] = sparseFormat(i,9*(I-1)+1:9*(I-1)+3,InvMahalanobis(Jri,cov));
                    [I_R2,J_R2,V_R2] = sparseFormat(i,9*(I-1)+7:9*(I-1)+9,InvMahalanobis(Jpi,cov));
                    [I_R3,J_R3,V_R3] = sparseFormat(i,9*(I)+1:9*(I)+3,InvMahalanobis(Jrip1,cov));
                    [I_R4,J_R4,V_R4] = sparseFormat(i,9*(I)+7:9*(I)+9,InvMahalanobis(Jpip1,cov));
                    [I_R5,J_R5,V_R5] = sparseFormat(i,mul*length(obj.states)+var_idx1,InvMahalanobis(Jlij,cov));
                    [I_R6,J_R6,V_R6] = sparseFormat(i,mul*length(obj.states)+var_idx2,InvMahalanobis(Jlijp1,cov));
                    [I_R7,J_R7,V_R7] = sparseFormat(i,mul*length(obj.states)+var_idx3,InvMahalanobis(Jlip1j,cov));
                    I_JacR(15*(i-1)+1:15*i) = [I_R1,I_R2,I_R3,I_R4,I_R5,I_R6,I_R7];
                    J_JacR(15*(i-1)+1:15*i) = [J_R1,J_R2,J_R3,J_R4,J_R5,J_R6,J_R7];
                    V_JacR(15*(i-1)+1:15*i) = [V_R1,V_R2,V_R3,V_R4,V_R5,V_R6,V_R7];
                end

                LM_jacR = sparse(I_JacR,J_JacR,V_JacR,blk_heightR,blk_width);

                LM_res = vertcat(LM_resL,LM_resR);
                LM_jac = vertcat(LM_jacL,LM_jacR);
            
            end
        end

        %% Retraction of Variables in Manifold
        function obj = retract(obj,delta,mode)
            n = length(obj.states);
            state_delta = delta(1:9*n);
            bias_ddelta = delta(9*n+1:15*n); % delta of delta value!

            obj.retractStates(state_delta);
            obj.retractBias(bias_ddelta,mode);
            
            if strcmp(obj.mode,'partial1')
                wsf_delta = delta(15*n+1:16*n);
                param_delta = delta(16*n+1:end);
                obj.retractWSF(wsf_delta,param_delta);
            elseif strcmp(obj.mode,'partial2')
                lane_delta = delta(15*n+1:end);
                obj.retractLane(lane_delta);
            elseif strcmp(obj.mode,'full')
                n_L = length(obj.I_L);
                n_R = length(obj.I_R);
                wsf_delta = delta(15*n+1:16*n);
                lane_delta = delta(16*n+1:16*n+n_L+n_R);
                param_delta = delta(16*n+n_L+n_R+1:end);
                obj.retractWSF(wsf_delta,param_delta);                
                obj.retractLane(lane_delta);
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
        function obj = retractWSF(obj,wsf_delta,param_delta)
            n = length(obj.states);
            for i=1:n                             
                obj.states{i}.WSF = obj.states{i}.WSF + wsf_delta(i);
            end           
            obj.params = obj.params + param_delta;
        end
        
        %% Lane Retraction
        function obj = retractLane(obj,lane_delta)
            % Add delta values to lane variables
            % Change format if optimizing 3D lane measurements
            n_L = length(obj.I_L); n_R = length(obj.I_R);
            deltaL = lane_delta(1:n_L);
            deltaR = lane_delta(n_L+1:end);
             
            for i=1:n_L
                I = obj.I_L(i); J = obj.J_L(i);
                obj.states{I}.left(2,J) = obj.states{I}.left(2,J) + deltaL(i);
            end

            for i=1:n_R
                I = obj.I_R(i); J = obj.J_R(i);
                obj.states{I}.right(2,J) = obj.states{I}.right(2,J) + deltaR(i);
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
            figure(2); spy(obj.opt.info); 
             
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
            
            disp('Extracting 2D Lane Point Covariance via Multivariate-Delta Method')
            for i=1:length(obj.states)
                R_idx = 9*(i-1)+1:9*(i-1)+3;
                P_idx = 9*(i-1)+7:9*(i-1)+9;
                state_cov = [obj.opt.cov(R_idx,R_idx), obj.opt.cov(R_idx,P_idx);
                             obj.opt.cov(P_idx,R_idx), obj.opt.cov(P_idx,P_idx)];
                lane_idx = obj.lane.state_idxs(i);
                
                R = obj.states{i}.R;

                for j=1:obj.lane.prev_num
                    LeftL_cov = diag([0, obj.lane.lystd(lane_idx,j)^2, obj.lane.lzstd(lane_idx,j)^2]);
                    RightL_cov = diag([0, obj.lane.rystd(lane_idx,j)^2, obj.lane.rzstd(lane_idx,j)^2]);
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
        function [Jrv,Jbg,Js,JL] = getWSSJac(Ri,Vi,Si,wi,bgi,bgdi,L)
            Jrv = zeros(3,6); % Horizontally Augmented for R and V
            Jrv(:,1:3) = 1/Si * skew(Ri' * Vi);
            Jrv(:,4:6) = 1/Si * Ri';
        
            Jbg = -1/Si * skew(L);
            Js = -1/Si^2 * (Ri' * Vi - skew(wi - bgi - bgdi) * L);
            JL = -1/Si * skew(wi - bgi - bgdi);
        end
        
        %% Lane Jacobians
        function [res_fixed,Jri,Jpi,Jrip1,Jpip1,Jlij,Jlijp1,Jlip1j] = getLaneJac(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J)
            
            res_fixed = getLaneRes(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J);
            
            delta = 1e-4; % Perturbation step
            % Ri Perturbation
            Jri = zeros(1,3);
            for i=1:3
                perturb_vec = zeros(3,1);
                perturb_vec(i) = delta;

                Ri_perturbedP = Ri * Exp_map(perturb_vec);
                Ri_perturbedM = Ri * Exp_map(-perturb_vec);
                res_perturbedP = getLaneRes(Ri_perturbedP,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J);
                res_perturbedM = getLaneRes(Ri_perturbedM,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J);
                Jri(i) = (res_perturbedP - res_perturbedM)/(2*delta);
            end
            
            Jpi = zeros(1,3);
            for i=1:3
                perturb_vec = zeros(3,1);
                perturb_vec(i) = delta;

                Pi_perturbedP = Pi + Ri * perturb_vec;
                Pi_perturbedM = Pi - Ri * perturb_vec;
                res_perturbedP = getLaneRes(Ri,Rip1,Pi_perturbedP,Pip1,Lij,Lijp1,Lip1j,J);
                res_perturbedM = getLaneRes(Ri,Rip1,Pi_perturbedM,Pip1,Lij,Lijp1,Lip1j,J);
                Jpi(i) = (res_perturbedP - res_perturbedM)/(2*delta);
            end

            Jrip1 = zeros(1,3);
            for i=1:3
                perturb_vec = zeros(3,1);
                perturb_vec(i) = delta;

                Rip1_perturbedP = Rip1 * Exp_map(perturb_vec);
                Rip1_perturbedM = Rip1 * Exp_map(-perturb_vec);
                res_perturbedP = getLaneRes(Ri,Rip1_perturbedP,Pi,Pip1,Lij,Lijp1,Lip1j,J);
                res_perturbedM = getLaneRes(Ri,Rip1_perturbedM,Pi,Pip1,Lij,Lijp1,Lip1j,J);
                Jrip1(i) = (res_perturbedP - res_perturbedM)/(2*delta);
            end
            
            Jpip1 = zeros(1,3);
            for i=1:3
                perturb_vec = zeros(3,1);
                perturb_vec(i) = delta;

                Pip1_perturbedP = Pip1 + Rip1 * perturb_vec;
                Pip1_perturbedM = Pip1 - Rip1 * perturb_vec;
                res_perturbedP = getLaneRes(Ri,Rip1,Pi,Pip1_perturbedP,Lij,Lijp1,Lip1j,J);
                res_perturbedM = getLaneRes(Ri,Rip1,Pi,Pip1_perturbedM,Lij,Lijp1,Lip1j,J);
                Jpip1(i) = (res_perturbedP - res_perturbedM)/(2*delta);
            end

            perturb_vec = delta;
            Lij_perturbedP = Lij + [0;perturb_vec;0];
            Lij_perturbedM = Lij - [0;perturb_vec;0];
            res_perturbedP = getLaneRes(Ri,Rip1,Pi,Pip1,Lij_perturbedP,Lijp1,Lip1j,J);
            res_perturbedM = getLaneRes(Ri,Rip1,Pi,Pip1,Lij_perturbedM,Lijp1,Lip1j,J);
            Jlij = (res_perturbedP - res_perturbedM)/(2*delta);

            perturb_vec = delta;
            Lijp1_perturbedP = Lijp1 + [0;perturb_vec;0];
            Lijp1_perturbedM = Lijp1 - [0;perturb_vec;0];
            res_perturbedP = getLaneRes(Ri,Rip1,Pi,Pip1,Lij,Lijp1_perturbedP,Lip1j,J);
            res_perturbedM = getLaneRes(Ri,Rip1,Pi,Pip1,Lij,Lijp1_perturbedM,Lip1j,J);
            Jlijp1 = (res_perturbedP - res_perturbedM)/(2*delta);

            perturb_vec = delta;
            Lip1j_perturbedP = Lip1j + [0;perturb_vec;0];
            Lip1j_perturbedM = Lip1j - [0;perturb_vec;0];
            res_perturbedP = getLaneRes(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j_perturbedP,J);
            res_perturbedM = getLaneRes(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j_perturbedM,J);
            Jlip1j = (res_perturbedP - res_perturbedM)/(2*delta);

%             Jri = (-Ri' * (skew(Lijp1) + lambda * Ri * skew(Lij - Lijp1)));
%             Jpi = Rip1' * Ri;
%             Jrip1 = (skew(Rip1' * (Pi + Ri * Lijp1 - Pip1)) + lambda * skew(Rip1' * Ri * (Lij - Lijp1)));
%             Jpip1 = -eye(3);
%             Jlam = Rip1' * Ri * (Lij - Lijp1);
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

        %% Data Association
        function idx = match2D(PC,Point)
            % Performs 2D index match of input Point Cloud and Point
            % If finds match, lower bound is returned.
            % Else, 0 is returned
            D = sqrt((PC(1,:) - Point(1)).^2 + (PC(2,:) - Point(2)).^2);
            
            [~,idx1] = min(D);
            if idx1 == 1
                % Test with idx = 1, 2
                d1 = D(idx1); d2 = D(idx1+1);
                d12 = sqrt((PC(:,idx1) - PC(:,idx1+1))' * (PC(:,idx1) - PC(:,idx1+1)));
                [~,max_idx] = max([d12,d1,d2]);
                if max_idx == 1
                    idx = 1;
                else
                    idx = 0;
                end
            elseif idx1 == size(PC,2)
                d1 = D(idx1); d2 = D(idx1-1);
                d12 = sqrt((PC(:,idx1) - PC(:,idx1-1))' * (PC(:,idx1) - PC(:,idx1-1)));
                [~,max_idx] = max([d12,d1,d2]);
                if max_idx == 1
                    idx = size(PC,2)-1;
                else
                    idx = 0;
                end
            else
                % Test with adjacent points
                dp = D(idx1-1); dn = D(idx1+1); dc = D(idx1);
                dpc = sqrt((PC(:,idx1) - PC(:,idx1-1))' * (PC(:,idx1) - PC(:,idx1-1)));
                dcn = sqrt((PC(:,idx1) - PC(:,idx1+1))' * (PC(:,idx1) - PC(:,idx1+1)));
                % Test previous
                [~,max_idxpc] = max([dpc,dp,dc]);
                [~,max_idxcn] = max([dcn,dc,dn]);

                if max_idxpc == 1 && max_idxcn == 1
                    if dp > dn
                        idx = idx1;
                    else
                        idx = idx1 - 1;
                    end
                else
                    if max_idxpc == 1
                        idx = idx1 - 1;
                    else
                        idx = idx1;
                    end
                end
                
            end

            % If there are multiple matches, pick the 'best' match with
            % smallest matching error.
            %
            % Perhaps change matching error to normalized error in the
            % future? (Reliability based) 
            % --> May be effective for inaccurate measurements
            
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
        
        %% Compute Error with different data length
        function E = computeError(p_ref,p_est)
            org_dlen = size(p_ref,2);
            dlen = size(p_est,2);
            E = zeros(1,dlen);

            ratio = (org_dlen - 1)/(dlen - 1);
            for i=1:dlen
                P_ref = p_ref(:,1+round((i-1)*ratio));
                P_est = p_est(:,i);
                E(i) = norm(P_ref - P_est);
            end
        end
        
    end
end

%% Lane Residuals
function res = getLaneRes(Ri,Rip1,Pi,Pip1,Lij,Lijp1,Lip1j,J)
%     rpy_i = dcm2rpy(Ri); psi_i = rpy_i(3);
%     rpy_ip1 = dcm2rpy(Rip1); psi_ip1 = rpy_ip1(3);
%     R2di = [cos(psi_i), -sin(psi_i);sin(psi_i), cos(psi_i)];
%     R2dip1 = [cos(psi_ip1), -sin(psi_ip1);sin(psi_ip1), cos(psi_ip1)];
%     
%     P2di = Pi(1:2); P2dip1 = Pip1(1:2);
%     L2dij = Lij(1:2); L2dijp1 = Lijp1(1:2); 
%     L2dip1j = Lip1j(2); % 1d variable

    Qijb = Rip1' * (Pi + Ri * Lij - (Rip1 * [10*(J-1);0;0] + Pip1));
    Qijp1b = Rip1' * (Pi + Ri * Lijp1 - (Rip1 * [10*(J-1);0;0] + Pip1));

    x1 = Qijb(1); y1 = Qijb(2);
    x2 = Qijp1b(1); y2 = Qijp1b(2);

    res = (x2*y1-x1*y2)/(x2-x1) - Lip1j(2);
end
