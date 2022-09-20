classdef optimizer < handle
% OPTIMIZER - Sensor Fusion Optimizer Module
% All input data are from the dataprocessor module
% [Usage Format]
% sol = OPTIMIZER(imu,gnss,lane,can,imu_bias,t,prior_covs'basic',options)
%
% [imu, gnss, lane, can, imu_bias, t] : Pre-processed data
% 
% [Mode] : 'basic', 'partial', 'full', '2-phase'
% * Basic: INS + GNSS Fusion
% * Partial: INS + GNSS + WSS Fusion
% * Full, 2-phase: INS + GNSS + WSS + Lane Detection Fusion (To be done)
% 
% [Options] : struct variable containing NLS optimization options
% * options.CostThres: Cost difference threshold
% * options.StepThres: Step size threshold
% * options.IterThres: Iteration limit
% * options.Algorithm: GN(Gauss-Newton),LM(Levenberg-Marquardt),TR(Trust-Region)
% ** GN and LM : no guarantee of convergence to local minima
% Recommended to use GN or TR 
% When using TR as solver algorithm, need to define parameters
% TR parameters at 'main.m' should work fine
% 
% [Methods] : using OPTIMIZER 
% * optimize() : Find optimized solution to SNLS problem
% * visualize() : Plot optimization results
% 
% Implemented by JinHwan Jeon, 2022

    properties
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
        opt = struct() % Optimized results
    end

    methods
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
            if strcmp(obj.mode,'partial')
                obj.params =  [1.5; 0; 1];  
            end

            % Perform Pre-integration
            obj.integrate(1:n-1);

            % Perform INS propagation
            obj.ins();

            % If mode is full, 2-phase, create initial value for lane
            % variables
        end
        
        %% Pre-integration using IMU Measurements, given IMU idxs
        function obj = integrate(obj,idxs)
            m = length(idxs);
            nbg_cov = obj.covs.imu.GyroscopeNoise;
            nba_cov = obj.covs.imu.AccelerometerNoise;
            n_cov = blkdiag(nbg_cov,nba_cov);

            for i=1:m
%                 disp(i)
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
                    [Ak, Bk] = getCoeff(delRik,delRkkp1,a_k,w_k,dt_k);

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
            if strcmp(obj.mode,'partial')
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
                state_.left = [0:10:10*(prev_num-1);
                               obj.lane.ly(lane_idx,1:prev_num);
                               obj.lane.lz(lane_idx,1:prev_num)];

                state_.right = [0:10:10*(prev_num-1);
                                obj.lane.ry(lane_idx,1:prev_num);
                                obj.lane.rz(lane_idx,1:prev_num)];

                obj.states = [obj.states {state_}];
            end

        end

        %% Optimize
        function obj = optimize(obj)
            
            % Optimize Sparse Nonlinear Least Squares in Batch Manner                     

            disp('[Optimization Starts...]')

            n = length(obj.imu)+1;

            if strcmp(obj.mode,'basic')
                obj.opt.x0 = zeros(15*n,1);
            elseif strcmp(obj.mode,'partial')
                obj.opt.x0 = zeros(16*n,1);
            end
            
            % Run optimization depending on algorithm options
            if strcmp(obj.opt.options.Algorithm,'GN')
                obj.GaussNewton();
            elseif strcmp(obj.opt.options.Algorithm,'LM')
                obj.LevenbergMarquardt();
            elseif strcmp(obj.opt.options.Algorithm,'TR')
                obj.TrustRegion();
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
                    osc_flag = DetectOsc(cost_stack);
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

            disp('[SNLS solver: Approximate Trust Region Method]')
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
                    string = 'Rejected';
                    flag = false;                    
                end
                
                step_size = norm(x0);
                
                str = sprintf(formatstr,i,cost,step_size,lambda,string);
                disp(str)
                
                lambda = UpdateLambda(lambda,Lu,Ld,flag);
                
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

            disp('[SNLS solver: Approximate Trust Region Method]')
            fprintf(' Iteration      f(x)        step        TR_radius    Acceptance\n');
            formatstr = ' %5.0f        %-10.3g   %-10.3g   %-10.3g    %s';
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
                h_gn = A \ b; % Gauss-Newton Step

                alpha = (b' * b)/(b' * A * b);
                h_gd = alpha * b; % Gradient Descent Step
                
                x0 = ComputeDogLeg(h_gn,h_gd,tr_rad);
                
                dummy_states = obj.states; % need to check if change in dummy state changes obj.states
                dummy_bias = obj.bias;
                dummy_imu = obj.imu;

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
                    string = 'Rejected';
                    flag = false;
                end
                
                step_size = norm(x0);
                
                str = sprintf(formatstr,i,cost,step_size,tr_rad,string);
                disp(str)
                
                tr_rad = UpdateDelta(rho,tr_rad,eta1,eta2,gamma1,gamma2);
                
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

        %% Optimization Cost Function
        function [res,jac] = cost_func(obj,x0)
            obj.retract(x0,'normal');
            [Pr_res,Pr_jac] = obj.CreatePrBlock();
            [MM_res,MM_jac] = obj.CreateMMBlock();
            [GNSS_res,GNSS_jac] = obj.CreateGNSSBlock();
            [WSS_res,WSS_jac] = obj.CreateWSSBlock();

            res = vertcat(Pr_res,MM_res,GNSS_res,WSS_res);
            jac = vertcat(Pr_jac,MM_jac,GNSS_jac,WSS_jac);            
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
            
            WSS_res = []; Params_res = [];
            if strcmp(obj.mode,'partial')
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
            
            Pr_res = [Veh_res; Bias_res; WSS_res; Params_res];
            Pr_jac = sparse(I,J,V,blk_height,blk_width);
        end
        
        %% IMU Residual and Jacobian
        function [MM_res,MM_jac] = CreateMMBlock(obj)
            n = length(obj.states);
            blk_width = 15*n;
%             blk_height = 15 * (n-1);
            grav = [0;0;-9.81];

            if strcmp(obj.mode,'partial')
                blk_width = blk_width + n;
                I_s = zeros(1,2*(n-1)); J_s = I_s; V_s = I_s;
                s_cov = obj.covs.imu.ScaleFactorNoise;
                WSF_res = zeros(n-1,1);

            elseif strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                % Add prev_num states
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

                [Ji,Jj,Jb] = getIMUJac(resR,Ri,Rj,Vi,Vj,Pi,Pj,bgdi,JdelRij_bg,...
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
                
                if strcmp(obj.mode,'partial')
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
            
            if strcmp(obj.mode,'partial')
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
                % Add prev_num states
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
                    % Add prev_num states
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
    
                    [Jrv,Jbg,Js] = getWSSJac(Ri,Vi,Si,wi,bgi,bgdi,L);

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

        %% Retraction
        function obj = retract(obj,delta,mode)
            
            n = length(obj.states);
            state_delta = delta(1:9*n);
            bias_ddelta = delta(9*n+1:15*n); % delta of delta value!

            obj.retractStates(state_delta);
            obj.retractBias(bias_ddelta,mode);
            
            if strcmp(obj.mode,'partial')
                wsf_delta = delta(15*n+1:16*n);
                obj.retractWSF(wsf_delta);
                
            elseif strcmp(obj.mode,'full') || strcmp(obj.mode,'2-phase')
                % Need to develop
                error("To be done")
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
                   bg_ = obj.bias{i}.bg + new_bgd;
                   bgd_ = zeros(3,1);
                   ba_ = obj.bias{i}.ba + new_bad;
                   bad_ = zeros(3,1);

                elseif strcmp(mode,'normal')
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
                plot(left_(1,:),left_(2,:));
                plot(right_(1,:),right_(2,:));
            end

            xlabel('Global X'); ylabel('Global Y');
            title('Optimized Trajectory 3D Comparison')
            legend([p_est,p_gnss],'Estimated Trajectory','GNSS Measurements')
            
            figure(2); 
%             snap_lla = [obj.snap.lat, obj.snap.lon, zeros(size(obj.snap.lat,1),1)];
            P_lla = enu2lla(P',obj.gnss.lla0,'ellipsoid');
            
            geoplot(P_lla(:,1),P_lla(:,2),'r.',...
                    obj.gnss.pos(:,1),obj.gnss.pos(:,2),'b.'); grid on; hold on;
%             geoplot(obj.snap.lat,obj.snap.lon,'c.','MarkerSize',8)

            geobasemap satellite

            title('Optimized Trajectory Comparison')
            legend('Estimated Trajectory','GNSS Measurements')

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
            legend('Vx_{b}','Vy_{b}','Vz_{b}','Wheel Speed')

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

            if strcmp(obj.mode,'partial')
                figure(8); 
                plot(obj.t,S); grid on;
                xlabel('Time(s)'); ylabel('Wheel Speed Scaling Factor')
                title('Optimized Wheel Speed Scaling Factor')
            end
        end

    end
end

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

function [Jrv,Jbg,Js] = getWSSJac(Ri,Vi,Si,wi,bgi,bgdi,L)
    Jrv = zeros(3,6); % Horizontally Augmented for R and V
    Jrv(:,1:3) = 1/Si * skew(Ri' * Vi);
    Jrv(:,4:6) = 1/Si * Ri';

    Jbg = -1/Si * skew(L);
    Js = -1/Si^2 * (Ri' * Vi - skew(wi - bgi - bgdi) * L);
%     Jl = -1/Si * skew(wi - bgi - bgdi);
%     Jl = Jl(:,1); % Only take the longitudinal lever arm as the optimization variable
end

function h_dl = ComputeDogLeg(h_gn,h_gd,tr_rad)
% Compute Dog-Leg Step for Approximate Trust Region Method
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

function Delta = UpdateDelta(rho,tr_rad,eta1,eta2,gamma1,gamma2)
% Update trust-region radius according to gain ratio value (rho)
    if rho >= eta2
        Delta = gamma2 * tr_rad;
    elseif rho < eta1
        Delta = gamma1 * tr_rad;
    else
        Delta = tr_rad;
    end
end

function Lambda = UpdateLambda(lambda,Lu,Ld,flag)
% Update Levenberg-Marquardt lambda value according to gain ration value
    if flag
        Lambda = max([lambda/Ld,1e-7]);
    else
        Lambda = min([lambda*Lu,1e7]);
    end
end

function osc_flag = DetectOsc(cost_stack)
% Detect oscillation by computing maximum deviation from average (recent 5 costs)
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
