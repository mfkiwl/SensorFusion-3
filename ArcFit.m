classdef ArcFit < handle
    %ArcFit NLS based simple Arc Fitting given data points
    %
    %   Given data points and initial arc parameters, this solver module
    %   finds the optimal arc parameter set.
    %   This module performs optimization for only given number of segments
    %   * No segment division/addition
    %   
    %   [NLS Models]
    %   Model 1: Arc-Point Measurement Model
    %   Model 2: Arc Initial/Final Point Anchoring Model
    %
    %   * The optimization process is very sensitive to covariance values.
    %   * Choosing inappropriate covariance weights may lead to singularity
    %   in extreme scenarios (Sharp Turn, etc)
    %
    %   [Optimization Variables]
    %   x0, y0, tau0, kappa1, ..., kappaN, L1, ..., LN
    %   --> Only arc parameters are to be optimized
    %   
    %   * 0m previewed lane points and remaining lane points are separately
    %   considered inside the measurement model
    %   
    %   * Since the covariance of lane points are given in 2D format,
    %   current implementation of ArcFit is not valid for usage.
    % 
    %   * ArcFit is deprecated
    %
    %   Implemented by JinHwan Jeon, 2022

    properties(Access = public)
        id
        params        
        points0
        points
        covs0
        covs        
        valid = false
        assoc0
        assoc
        arc_centers
        arc_nodes
        precomp_jac = struct()
        validity
%         options = struct()
        opt = struct()
    end

    methods(Access = public)        
        %% Constructor
        function obj = ArcFit(params,points0,points,covs0,covs,id)
            obj.id = id;
            obj.params = params;
            obj.points0 = points0;
            obj.points = points;
            obj.covs0 = covs0;
            obj.covs = covs;
            obj.assoc = zeros(1,size(obj.points,2));
            
%             obj.options.CostThres = 1e-1;
%             obj.options.StepThres = 1e-10;
%             obj.options.IterThres = inf;
%             obj.options.TR = struct();            
%             obj.options.TR.eta1 = 0.6;
%             obj.options.TR.eta2 = 0.9;
%             obj.options.TR.gamma1 = 0.1;
%             obj.options.TR.gamma2 = 2;
%             obj.options.TR.thres = 1e-10; % Trust Region Radius Threshold
%             obj.options.minL = 5; % Minimum arc length constraint
        end
            
        %% Retrieve params 
        function params = getParams(obj)
            params = obj.params;
        end
        
        %% Optimize
        function obj = optimize(obj)
            obj.associate0();
            obj.associate();
            disp(['[Performing Optimization for Segment ID: ',num2str(obj.id),']']) 
            X0 = [obj.params.x0; obj.params.y0; obj.params.tau0; obj.params.kappa'; obj.params.L'];
%             obj.TrustRegion();
            n = length(obj.params.kappa);
            lb = [repmat(-inf,1,3+n),zeros(1,n)];
            ub = [];
            options = optimoptions('lsqnonlin', ...
                                   'UseParallel',true, ...
                                   'Display','iter-detailed', ...
                                   'MaxFunctionEvaluations',inf, ...    
                                   'MaxIterations',inf, ...
                                   'FiniteDifferenceType','central');
            X = lsqnonlin(@obj.cost_func,X0,lb,ub,options);
            
            obj.params.x0 = X(1);
            obj.params.y0 = X(2);
            obj.params.tau0 = X(3);
            obj.params.kappa = X(3+1:3+n)';
            obj.params.L = X(3+n+1:end)';
        end
        
        %% Visualize: One Segment Optimization
        function visualize(obj)
            figure(1);
            p_est = plot(obj.points(1,:),obj.points(2,:),'r.'); hold on; grid on; axis equal;
            n = length(obj.params.kappa);
            heading = obj.params.tau0;
            SegPoints = [obj.params.x0;
                         obj.params.y0];
            for i=1:n
                kappa = obj.params.kappa(i); L = obj.params.L(i);
                headingPrev = heading;
                heading = heading + kappa * L;
                headingCurr = heading;

                heading_ = linspace(headingPrev,headingCurr,1e3);
                addedSegPoints = SegPoints(:,end) + 1/kappa * [sin(heading_) - sin(headingPrev);
                                                               -cos(heading_) + cos(headingPrev)];
                SegPoints = [SegPoints addedSegPoints];
                
                p_node = plot(addedSegPoints(1,1),addedSegPoints(2,1),'co');
                plot(addedSegPoints(1,end),addedSegPoints(2,end),'co');
            end
            p_lane = plot(SegPoints(1,:),SegPoints(2,:),'k-');            
            xlabel('Global X(m)'); ylabel('Global Y(m)');
            title('Optimized Vehicle Trajectory and Lane Segments');
            legend([p_est,p_lane,p_node], ...
                   'Input Lane Points', ...
                   'Optimized Arc Spline', ...
                   'Lane Sub-segment');
        end
        
    end
    methods(Access = private)  
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
            n = length(obj.params.kappa);
            x0 = zeros(2*n+3,1);
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
                obj.opt.flags = [obj.opt.flags abs(prev_cost - cost) > obj.options.CostThres];
                obj.opt.flags = [obj.opt.flags step_size > obj.options.StepThres];
                obj.opt.flags = [obj.opt.flags i < obj.options.IterThres];
                
                % Check for oscillation around the local minima
                if length(cost_stack) >= 5
                    osc_flag = obj.DetectOsc(cost_stack);
                else
                    osc_flag = false;
                end

                if length(find(obj.opt.flags)) ~= length(obj.opt.flags) % If any of the criterion is not met, end loop
                    
                    obj.retract(x0);

                    disp('[Optimization Finished...]')
                    idx = find(~obj.opt.flags,1);
                    
                    if idx == 1
                        disp(['Current cost difference ',num2str(abs(prev_cost-cost)),' is below threshold: ',num2str(obj.options.CostThres)])
                    elseif idx == 2
                        disp(['Current step size ',num2str(step_size),' is below threshold: ',num2str(obj.options.StepThres)])
                    elseif idx == 3
                        disp(['Current iteration number ',num2str(i),' is above threshold: ',num2str(obj.options.IterThres)])
                    end

                    break;
                elseif osc_flag
                    obj.retract(x0);

                    disp('[Optimization Finished...]')
                    disp('Oscillation about the local minima detected')
                    break;
                else
                    i = i + 1;
                    prev_cost = cost;
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
            fprintf(' Iteration      f(x)          step        TR_radius    Acceptance\n');
            formatstr = ' %5.0f        %-12.6g   %-10.3g   %-10.3g    %s';
            
            n = length(obj.params.kappa);
            x0 = zeros(2*n+3,1);

            tr_rad = 10; % Initial Trust Region radius
            [res,jac] = obj.cost_func(x0);
            prev_cost = res' * res;
            str = sprintf(formatstr,0,prev_cost,norm(x0),tr_rad,'Init');
            disp(str)

            i = 1;

            eta1 = obj.options.TR.eta1;
            eta2 = obj.options.TR.eta2;
            gamma1 = obj.options.TR.gamma1;
            gamma2 = obj.options.TR.gamma2;

            while true
                A = jac' * jac; b = -jac' * res; 
                
                h_gn = A \ b; % Gauss-Newton Step

                alpha = (b' * b)/(b' * A * b);
                h_gd = alpha * b; % Gradient Descent Step                

                x0 = obj.ComputeDogLeg(h_gn,h_gd,tr_rad);                
                
                dummy_params = obj.params;
                dummy_assoc = obj.assoc;
                dummy_jac = obj.precomp_jac;

                if ~isempty(find(isnan(h_gn),1)) || ~isempty(find(isnan(h_gd),1))
                    obj.opt.jac = jac;
                    obj.opt.info = A;                    

                    error('Singular Matrix')
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
                    % Current Step Rejected, using dummy variables
                    obj.params = dummy_params;   
                    obj.assoc = dummy_assoc;
                    obj.precomp_jac = dummy_jac;
                    string = 'Rejected';
                    flag = false;
                end
                
                step_size = norm(x0);
                
                str = sprintf(formatstr,i,cost,step_size,tr_rad,string);
                disp(str)
                
                tr_rad = obj.UpdateDelta(rho,tr_rad,eta1,eta2,gamma1,gamma2);
                
                if flag % Check Ending Criterion for Accepted Steps
                    flags = [];
                    flags = [flags abs(ared) > obj.options.CostThres];
                    flags = [flags step_size > obj.options.StepThres];
                    flags = [flags i < obj.options.IterThres];
                    flags = [flags tr_rad > obj.options.TR.thres];

                    if length(find(flags)) ~= length(flags) % If any of the criterion is not met, end loop
                        
                        obj.retract(x0);
    
                        disp('[Optimization Finished...]')
                        idx = find(~flags,1);
                        
                        if idx == 1
                            disp(['Current cost difference ',num2str(abs(ared)),' is below threshold: ',num2str(obj.options.CostThres)])
                        elseif idx == 2
                            disp(['Current step size ',num2str(step_size),' is below threshold: ',num2str(obj.options.StepThres)])
                        elseif idx == 3
                            disp(['Current iteration number ',num2str(i),' is above threshold: ',num2str(obj.options.IterThres)])
                        elseif idx == 4
                            disp(['Current trust region radius ',num2str(tr_rad),' is below threshold: ',num2str(obj.options.TR.thres)])
                        end
    
                        break;
                    end

                    i = i + 1;
                else
                    if tr_rad < obj.options.TR.thres
                        disp(['Current trust region radius ',num2str(tr_rad),' is below threshold: ',num2str(obj.options.TR.thres)])
                        break;
                    end
                end
            end
        end
                
        %% NLS Cost Function Evaluation
        function res = cost_func(obj,x0)
%             obj.retract(x0);
            n = length(obj.params.kappa);
            initParams = x0(1:3);
            kappa = x0(3+1:3+n);
            L = x0(3+n+1:end);
            ME_res = obj.CreateMEBlock(initParams,kappa,L);
            AM_res = obj.CreateAMBlock(initParams,kappa,L);
%             [AL_res,AL_jac] = obj.CreateALBlock();
            res = vertcat(ME_res,AM_res);
%             jac = vertcat(ME_jac,AM_jac);
        end

        %% Create Measurement Block
        function res = CreateMEBlock(obj,initParams,kappa,L)

            n = length(obj.params.kappa);
%             blk_width = 2*n+3;
            % 0m Previewed Measurements
            blk_height0 = size(obj.points0,2);            
            res0 = zeros(blk_height0,1);
            
%             I0 = zeros(1,3*size(obj.points0,2)+2*sum(obj.assoc0));
%             J0 = I0; V0 = I0;

%             jac_cnt0 = 0;
            for i=1:size(obj.points0,2)
                SegIdx = obj.assoc0(i);
                Point = obj.points0(:,i);
%                 cov = obj.covs0(i);
                cov = 0.1^2;
%                 [init_jac,kappa_jac,L_jac,anchored_res] = obj.MEjac(SegIdx,Point,cov);
                
                res0(i) = InvMahalanobis(obj.MEres(initParams,kappa,L,SegIdx,Point),cov);
%                 [I_P,J_P,V_P] = sparseFormat(i,1:3,init_jac);
%                 [I_k,J_k,V_k] = sparseFormat(i,3+1:3+SegIdx,kappa_jac);
%                 [I_l,J_l,V_l] = sparseFormat(i,3+n+1:3+n+SegIdx,L_jac);
%                 I0(jac_cnt0+1:jac_cnt0+3+2*SegIdx) = [I_P, I_k, I_l];
%                 J0(jac_cnt0+1:jac_cnt0+3+2*SegIdx) = [J_P, J_k, J_l];
%                 V0(jac_cnt0+1:jac_cnt0+3+2*SegIdx) = [V_P, V_k, V_l];
%                 jac_cnt0 = jac_cnt0 + 3 + 2 * SegIdx;
            end
%             jac0 = sparse(I0,J0,V0,blk_height0,blk_width);

            % Remaining Previewed Measurements
            blk_heightR = nnz(obj.assoc);
            resR = zeros(blk_heightR,1);

%             IR = zeros(1,3*blk_heightR + 2*sum(obj.assoc));
%             JR = IR; VR = IR;

%             jac_cntR = 0;
            cntR = 1;
            for i=1:size(obj.points,2)
                SegIdx = obj.assoc(i);
                if SegIdx > 0
                    Point = obj.points(:,i);
%                     cov = obj.covs(i);
                    cov = 0.3^2;
%                     [init_jac,kappa_jac,L_jac,anchored_res] = obj.MEjac(SegIdx,Point,cov);
                    
                    resR(cntR) = InvMahalanobis(obj.MEres(initParams,kappa,L,SegIdx,Point),cov);
%                     [I_P,J_P,V_P] = sparseFormat(cntR,1:3,init_jac);
%                     [I_k,J_k,V_k] = sparseFormat(cntR,3+1:3+SegIdx,kappa_jac);
%                     [I_l,J_l,V_l] = sparseFormat(cntR,3+n+1:3+n+SegIdx,L_jac);
%                     IR(jac_cntR+1:jac_cntR+3+2*SegIdx) = [I_P, I_k, I_l];
%                     JR(jac_cntR+1:jac_cntR+3+2*SegIdx) = [J_P, J_k, J_l];
%                     VR(jac_cntR+1:jac_cntR+3+2*SegIdx) = [V_P, V_k, V_l];
%                     jac_cntR = jac_cntR + 3 + 2 * SegIdx;
                    cntR = cntR + 1;                    
                end
            end
%             jacR = sparse(IR,JR,VR,blk_heightR,blk_width);
%             resR = []; jacR = [];
            res = vertcat(res0,resR);   
%             jac = vertcat(jac0,jacR);
        end
        
        %% Create Measurement Jacobian 
        function [init_jac,kappa_jac,L_jac,anchored_res] = MEjac(obj,SegIdx,Point,cov)
            
            kappa = obj.params.kappa(1:SegIdx); 
            
            Xc = obj.arc_centers(:,SegIdx);
            anchored_res = sqrt((Xc - Point)' * (Xc - Point)) - abs(1/kappa(SegIdx));            
            init_jac = zeros(1,3);
            kappa_jac = zeros(1,SegIdx);
            L_jac = zeros(1,SegIdx);

%             % Algebraic Representation for fast jacobian computation
%             % Pre-computed Jacobian values are computed every iteration,
%             % after data association: Need to fix
%             xd = Point(1); yd = Point(2);
%             xc = obj.arc_centers(1,SegIdx); yc = obj.arc_centers(2,SegIdx);
%             
%             delX = (xc - xd); delY = (yc - yd);
%             A = 1/sqrt(delX^2 + delY^2);
%             
%             n = length(kappa);
%             init_jac(1) = A * delX;
%             init_jac(2) = A * delY;
%             init_jac(3) = A * (delX * obj.precomp_jac.Xc(SegIdx,3) + delY * obj.precomp_jac.Yc(SegIdx,3));
%             
%             for i=1:SegIdx
%                 
%                 kappa_jac(i) = A * (delX * obj.precomp_jac.Xc(SegIdx,3+i) + delY * obj.precomp_jac.Yc(SegIdx,3+i));
%                 L_jac(i) = A * (delX * obj.precomp_jac.Xc(SegIdx,3+n+i) + delY * obj.precomp_jac.Yc(SegIdx,3+n+i));
%                 if i == SegIdx
%                     kappa_jac(i) = kappa_jac(i) + sign(kappa(SegIdx))/kappa(SegIdx)^2;
%                 end
%             end
% 
%             init_jac = InvMahalanobis(init_jac,cov);
%             kappa_jac = InvMahalanobis(kappa_jac,cov);
%             L_jac = InvMahalanobis(kappa_jac,cov);
%             anchored_res = InvMahalanobis(anchored_res,cov);
            
            % Numerical Jacobian Computation
            
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa(1:SegIdx); L = obj.params.L(1:SegIdx);
            eps = 1e-8; eps_kappa = 1e-10; 
            minL = obj.options.minL;
            % init_jac
            for i=1:3
                initParams_ptb_vec = zeros(1,3);
                initParams_ptb_vec(i) = eps;
                initParams_ptb = initParams + initParams_ptb_vec;
                res_ptb = obj.MEres(initParams_ptb,kappa,L,SegIdx,Point);
                init_jac(i) = (res_ptb - anchored_res)/eps;
            end

            % kappa_jac, L_jac
            for i=1:SegIdx
                kappa_ptb_vec = zeros(1,SegIdx);
                L_ptb_vec = zeros(1,SegIdx);
                kappa_ptb_vec(i) = eps_kappa;
                L_ptb_vec(i) = eps;
                kappa_ptb = kappa + kappa_ptb_vec;
                L_ptb = ((L - minL).^(1/2) + L_ptb_vec).^2 + minL;
                res_ptbK = obj.MEres(initParams,kappa_ptb,L,SegIdx,Point);
                res_ptbL = obj.MEres(initParams,kappa,L_ptb,SegIdx,Point);
                kappa_jac(i) = (res_ptbK - anchored_res)/eps_kappa;
                L_jac(i) = (res_ptbL - anchored_res)/eps;
            end           

            init_jac = InvMahalanobis(init_jac,cov);
            kappa_jac = InvMahalanobis(kappa_jac,cov);
            L_jac = InvMahalanobis(L_jac,cov);

            anchored_res = InvMahalanobis(anchored_res,cov);
        end
        
        %% Arc Measurement Residual
        function res = MEres(obj,initParams,kappa,L,SegIdx,Point)
            centerPoints = obj.propCenter(initParams,kappa,L);
            Xc = centerPoints(:,SegIdx);
            res = sqrt((Xc - Point)' * (Xc - Point)) - abs(1/kappa(SegIdx));
        end

        %% Create Anchor Measurement Block
        function res = CreateAMBlock(obj,initParams,kappa,L)
            blk_height = 2*(size(obj.params.bnds,1)+1);
%             blk_height = 4;
            n = length(obj.params.kappa);
%             blk_width = 2 * n + 3;
            res = zeros(blk_height,1);
%             jac = zeros(blk_height,blk_width);
            
            
            % First
%             [init_jac,~,~,anchored_res] = obj.AMjac(0);            
%             res(1:2) = InvMahalanobis(obj.AMres(initParams,kappa,L,0),cov);
%             jac(1:2,1:3) = init_jac; 
%             jac(1:2,3+1:3+n) = kappa_jac;
%             jac(1:2,3+n+1:end) = L_jac;

            % Last
%             [init_jac,kappa_jac,L_jac,anchored_res] = obj.AMjac(n);
%             res(3:4) = InvMahalanobis(obj.AMres(initParams,kappa,L,n),cov);
%             jac(3:4,1:3) = init_jac; 
%             jac(3:4,3+1:3+n) = kappa_jac;
%             jac(3:4,3+n+1:end) = L_jac;
%             jac = sparse(jac);
            
%             if SubSegIdx == 0 || SubSegIdx == length(kappa)
%                 cov = diag([1e-4, 1e-4]); 
%             else
%                 cov = diag([0.5^2, 0.5^2]);
%             end

            for i=0:n
%                 [init_jac,kappa_jac,L_jac,anchored_res] = obj.AMjac(i);
                if i == 0 || i == n
                    cov = diag([1e-4, 1e-4]); 
                else
                    cov = diag([0.5^2, 0.5^2]);
                end
                res(2*i+1:2*i+2) = InvMahalanobis(obj.AMres(initParams,kappa,L,i),cov);
%                 jac(2*i+1:2*i+2,1:3) = init_jac;
%                 jac(2*i+1:2*i+2,3+1:3+n) = kappa_jac;
%                 jac(2*i+1:2*i+2,3+n+1:end) = L_jac;
            end
%             jac = sparse(jac);
        end

        %% Create Anchor Measurement Jacobian
        function [init_jac,kappa_jac,L_jac,anchored_res] = AMjac(obj,SubSegIdx)
            eps = 1e-8; eps_kappa = 1e-10;
            minL = obj.options.minL;
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            
            kappa = obj.params.kappa; L = obj.params.L;
            anchored_res = obj.AMres(initParams,kappa,L,SubSegIdx);

            init_jac = zeros(2,3);
            kappa_jac = zeros(2,length(kappa));
            L_jac = zeros(2,length(L));

            if SubSegIdx == 0 || SubSegIdx == length(kappa)
                cov = diag([1e-4, 1e-4]); 
            else
                cov = diag([0.5^2, 0.5^2]);
            end
            % Numerical Jacobian Computation
            % 1: Perturbation of initial parameters 
            for i=1:3
                initParams_ptb_vec = zeros(1,3);
                initParams_ptb_vec(i) = eps;
                initParams_ptbP = initParams + initParams_ptb_vec;
                res_ptb = obj.AMres(initParams_ptbP,kappa,L,SubSegIdx);
                init_jac(:,i) = (res_ptb - anchored_res)/eps;
            end
            
            init_jac = InvMahalanobis(init_jac,cov);

            % 2: Perturbation of Segment parameters 
            for i=1:SubSegIdx
                kappa_ptb_vec = zeros(1,length(kappa));
                kappa_ptb_vec(i) = eps_kappa;
                kappa_ptb = kappa + kappa_ptb_vec;

                L_ptb_vec = zeros(1,length(L));
                L_ptb_vec(i) = eps;
                L_ptb = ((L - minL).^(1/2) + L_ptb_vec).^2 + minL;
%                 L_ptb = L + L_ptb_vec;

                res_ptbK = obj.AMres(initParams,kappa_ptb,L,SubSegIdx);
                res_ptbL = obj.AMres(initParams,kappa,L_ptb,SubSegIdx);
                kappa_jac(:,i) = (res_ptbK - anchored_res)/eps_kappa;
                L_jac(:,i) = (res_ptbL - anchored_res)/eps;
            end
            
            kappa_jac = InvMahalanobis(kappa_jac,cov);
            L_jac = InvMahalanobis(L_jac,cov);

            anchored_res = InvMahalanobis(anchored_res,cov);
            
            % Algebraic Jacobian Computation
%             if SubSegIdx == 0
%                 init_jac = [1, 0, 0;
%                             0, 1, 0];
%                 init_jac = InvMahalanobis(init_jac,cov);
%             else
%                 init_jac(1,:) = obj.precomp_jac.Xn(SubSegIdx,1:3);
%                 init_jac(2,:) = obj.precomp_jac.Yn(SubSegIdx,1:3);
%                 init_jac = InvMahalanobis(init_jac,cov);
% 
%                 kappa_jac(1,:) = obj.precomp_jac.Xn(SubSegIdx,3+1:3+length(kappa));
%                 kappa_jac(2,:) = obj.precomp_jac.Yn(SubSegIdx,3+1:3+length(kappa));
%                 kappa_jac = InvMahalanobis(kappa_jac,cov);
% 
%                 L_jac(1,:) = obj.precomp_jac.Xn(SubSegIdx,3+length(kappa)+1:end);
%                 L_jac(2,:) = obj.precomp_jac.Yn(SubSegIdx,3+length(kappa)+1:end);
%                 L_jac = InvMahalanobis(L_jac,cov);
%             end
% 
%             anchored_res = InvMahalanobis(anchored_res,cov);
        end

        %% Anchor Measurement Residual
        function res = AMres(obj,initParams,kappa,L,SubSegIdx)
            nodePoints = obj.propNode(initParams,kappa,L);
            if SubSegIdx == 0   
                point = obj.points0(:,1);
                X = nodePoints(:,1);
            else
                bnds = obj.params.bnds;
                idx = bnds(SubSegIdx,2);
                point = obj.points0(:,idx);         
                X = nodePoints(:,SubSegIdx+1);    
            end                        
            res = X - point;
        end
        
        %% Arc Length Minimization Block
        function [res,jac] = CreateALBlock(obj)
            res = sum(obj.params.L);
            n = length(obj.params.L);
            jac = zeros(1,2*n+3);
            for i=1:n
                jac(3+n+i) = 1;
            end
        end

        %% Retract delta values
        function obj = retract(obj,dX)
            n = (length(dX)-3)/2;
            minL = obj.options.minL;
            obj.params.x0 = obj.params.x0 + dX(1);
            obj.params.y0 = obj.params.y0 + dX(2);
            obj.params.tau0 = obj.params.tau0 + dX(3);
            obj.params.kappa = obj.params.kappa + dX(4:3+n)';
            obj.params.L = ((obj.params.L - minL).^(1/2) + dX(4+n:end)').^2 + minL;
%             obj.params.L = obj.params.L + dX(4+n:end)';
            
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa; L = obj.params.L;
            obj.arc_nodes = obj.propNode(initParams,kappa,L);
            obj.arc_centers = obj.propCenter(initParams,kappa,L);
%             Arc Node and Center points are calculated every data
%             association
%             No need to save as dummy variables for Trust Region Method            

            % Perform data association
%             obj.associate(); 

            % Pre Compute Measurement Jacobians (Chain Rule)
%             obj.precomputeMEJac();          

            % Pre Compute Anchoring Model Jacobians
%             obj.precomputeAMJac();
        end
        
        %% Data Association 0m
        function obj = associate0(obj)
            obj.assoc0 = zeros(1,size(obj.points0,2));
            % If state idx is used, then association is fixed throughout
            % the whole optimization process
            bnds = obj.params.bnds;
            for i=1:size(bnds,1)
                lb = bnds(i,1); ub = bnds(i,2);
                obj.assoc0(lb:ub) = i;
            end
        end

        %% Data Association
        function obj = associate(obj)
            %ASSOCIATE Data Association 
            % Matches which point belongs to which segment
            obj.assoc = zeros(1,size(obj.points,2));            
            
            % Use arc centers to divide points 
            
%             cpydPoints = obj.points(1:2,:);
%             idxCounter = 1;
%             for i=1:n
%                 if i == n
%                     C1 = obj.arc_centers(:,i-1);
%                     C2 = obj.arc_centers(:,i);
%                 else
%                     C1 = obj.arc_centers(:,i);
%                     C2 = obj.arc_centers(:,i+1);
%                 end                
%                 
%                 x_init = cpydPoints(1,1); y_init = cpydPoints(2,1);
%                 slope = (C2(2) - C1(2))/(C2(1) - C1(1));
%                 val = y_init - (slope * (x_init - C1(1)) + C1(2));
%                 new_val = val;
%                 
%                 j = 0;
% %                 disp(['Segment Idx: ',num2str(i)])
% %                 disp(['Init Func val: ',num2str(val)])
%                 while val * new_val >= 0
%                     j = j + 1;
%                     if j > size(cpydPoints,2)
%                         break;
%                     end
%                     x_point = cpydPoints(1,j); y_point = cpydPoints(2,j);
%                     new_val = y_point - (slope * (x_point - C1(1)) + C1(2));
% %                     disp(['Search Idx: ',num2str(idxCounter+j-1)])
%                 end
%                 obj.assoc(idxCounter:idxCounter+j-2) = i * ones(1,j-1);
%                 if j > size(cpydPoints,2)
%                     break;
%                 end
%                 cpydPoints = cpydPoints(:,j:end);
%                 idxCounter = idxCounter + j-1;               
%             end
            
%             error('1');
                
            for i=1:size(obj.points,2)
                Px = obj.points(1,i); Py = obj.points(2,i);
                P_refx = obj.points0(1,:); P_refy = obj.points0(2,:);
                D = sqrt((P_refx - Px).^2 + (P_refy - Py).^2);
                [~,idx] = min(D);
                obj.assoc(i) = obj.assoc0(idx);
            end

%             % Find Match
%             for i=1:length(obj.points)
%                 % If not matched, 0 is returned
%                 obj.assoc(i) = obj.findMatch(obj.arc_nodes,obj.points(:,i));
%             end

            % Apply zero-order hold for zeros
%             obj.assoc = obj.ZOH(obj.assoc);                       
        end
                
        %% Find match for data association
        function idx = findMatch(obj,nodePoints,point)
            n = size(nodePoints,2);
            cand_idx = [];
            N_dist = [];
            kappa = obj.params.kappa;            

            for i=1:n-1
                P1 = nodePoints(:,i); P2 = nodePoints(:,i+1);
                x1 = P1(1); y1 = P1(2); x2 = P2(1); y2 = P2(2);
                xp = point(1); yp = point(2);
                slope = (y2-y1)/(x2-x1);

                x_vert = (xp + x1 * slope^2 + (yp-y1)*slope)/(1 + slope^2);
                y_vert = slope * (x_vert - x1) + y1;
                X_vert = [x_vert; y_vert];

                Dp = (P1 - P2)' * (P1 - P2);
                D1 = (P1 - X_vert)' * (P1 - X_vert);
                D2 = (P2 - X_vert)' * (P2 - X_vert);

                [~,max_idx] = max([Dp,D1,D2]);               

                if max_idx == 1
                    cand_idx = [cand_idx i];  
                    Xc = obj.arc_centers(:,i);                    
                    N_dist = [N_dist sqrt((Xc - point)' * (Xc - point)) - abs(1/kappa(i))];
                end
            end

            % If there are many possible matches, compare each match
            % Using the vertical distance
            
            if isempty(cand_idx)
                idx = 0;
            else
                [~,min_idx] = min(N_dist);
                idx = cand_idx(min_idx);
            end
        end

        %% Precompute Measurement Jacobian for Chain Rule
        function obj = precomputeMEJac(obj)
            % Iteratively pre compute jacobian terms for fast jacobian
            % computation
            % Computes Jacobian of all parameters w.r.t. Xc, Yc variables
            % (Xc, Yc) are center of arc segments 
            % 
            n = length(obj.params.kappa);           
            obj.precomp_jac.Xc = zeros(n,2*n+3);
            obj.precomp_jac.Yc = zeros(n,2*n+3);
            
            kappa = obj.params.kappa;
            L = obj.params.L;
            heading = obj.params.tau0;            

            for i=1:n
                % I : Index of Matched Sub-segment index
                % J : Index of Sub-segment parameter of interest

                if i == 1
                    obj.precomp_jac.Xc(1,1) = 1;
                    obj.precomp_jac.Xc(1,3) = -1/kappa(1) * cos(heading);
                    obj.precomp_jac.Xc(1,4) = 1/kappa(1)^2 * sin(heading);

                    obj.precomp_jac.Yc(1,2) = 1;
                    obj.precomp_jac.Yc(1,3) = -1/kappa(1) * sin(heading);
                    obj.precomp_jac.Yc(1,4) = -1/kappa(1)^2 * cos(heading);
                else
                    obj.precomp_jac.Xc(i,1) = obj.precomp_jac.Xc(i-1,1);
                    obj.precomp_jac.Xc(i,3) = obj.precomp_jac.Xc(i-1,3) + (1/kappa(i-1) - 1/kappa(i)) * cos(heading);

                    obj.precomp_jac.Yc(i,2) = obj.precomp_jac.Yc(i-1,2);
                    obj.precomp_jac.Yc(i,3) = obj.precomp_jac.Yc(i-1,3) + (1/kappa(i-1) - 1/kappa(i)) * sin(heading);

                    for j=1:i
                        % Kappa
                        if j == i-1
                            obj.precomp_jac.Xc(i,3+i-1) = obj.precomp_jac.Xc(i-1,3+i-1) + L(i-1) * (1/kappa(i-1) - 1/kappa(i)) * cos(heading) - 1/kappa(i-1)^2 * sin(heading);
                            obj.precomp_jac.Yc(i,3+i-1) = obj.precomp_jac.Yc(i-1,3+i-1) + L(i-1) * (1/kappa(i-1) - 1/kappa(i)) * sin(heading) + 1/kappa(i-1)^2 * cos(heading);
                        elseif j == i
                            obj.precomp_jac.Xc(i,3+i) = 1/kappa(i)^2 * sin(heading);
                            obj.precomp_jac.Yc(i,3+i) = -1/kappa(i)^2 * cos(heading);
                        else
                            obj.precomp_jac.Xc(i,3+j) = obj.precomp_jac.Xc(i-1,3+j) + L(j) * (1/kappa(i-1) - 1/kappa(i)) * cos(heading);
                            obj.precomp_jac.Yc(i,3+j) = obj.precomp_jac.Yc(i-1,3+j) + L(j) * (1/kappa(i-1) - 1/kappa(i)) * sin(heading);
                        end
                        % L
                        if j ~= i
                            obj.precomp_jac.Xc(i,3+n+j) = obj.precomp_jac.Xc(i-1,3+n+j) + kappa(j) * (1/kappa(i-1) - 1/kappa(i)) * cos(heading);
                            obj.precomp_jac.Yc(i,3+n+j) = obj.precomp_jac.Yc(i-1,3+n+j) + kappa(j) * (1/kappa(i-1) - 1/kappa(i)) * sin(heading);
                        end
                    end                    
                end

                heading = heading + kappa(i) * L(i);
            end
        end
        
        %% Precompute Anchoring Model Jacobian
        function obj = precomputeAMJac(obj)
            heading = obj.params.tau0;
            kappa = obj.params.kappa; L = obj.params.L;
            n = length(kappa);
            obj.precomp_jac.Xn = zeros(n,2*n+3);
            obj.precomp_jac.Yn = zeros(n,2*n+3);
            
            for i=1:n
                if i == 1 % x1, y1 : 1 step propagated point
                    obj.precomp_jac.Xn(1,1) = 1;
                    obj.precomp_jac.Xn(1,3) = 1/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));
                    obj.precomp_jac.Xn(1,4) = L(i)/kappa(i) * cos(heading + kappa(i) * L(i)) - 1/kappa(i)^2 * (sin(heading + kappa(i) * L(i)) - sin(heading));
                    obj.precomp_jac.Xn(1,3+n+1) = cos(heading + kappa(i) * L(i));

                    obj.precomp_jac.Yn(1,2) = 1;
                    obj.precomp_jac.Yn(1,3) = 1/kappa(i) * (sin(heading + kappa(i) * L(i)) - sin(heading));
                    obj.precomp_jac.Yn(1,4) = L(i)/kappa(i) * sin(heading + kappa(i) * L(i)) + 1/kappa(i)^2 * (cos(heading + kappa(i) * L(i)) - cos(heading));
                    obj.precomp_jac.Yn(1,3+n+1) = sin(heading + kappa(i) * L(i));
                else
                    obj.precomp_jac.Xn(i,1) = obj.precomp_jac.Xn(i-1,1);
                    obj.precomp_jac.Xn(i,3) = obj.precomp_jac.Xn(i-1,1) + 1/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));

                    obj.precomp_jac.Yn(i,2) = obj.precomp_jac.Yn(i-1,2);
                    obj.precomp_jac.Yn(i,3) = obj.precomp_jac.Yn(i-1,3) + 1/kappa(i) * (sin(heading + kappa(i) * L(i)) - sin(heading));
                    
                    for j=1:i
                        if j ~= i
                            obj.precomp_jac.Xn(i,3+j) = obj.precomp_jac.Xn(i-1,3+j) + 1/kappa(i) * L(j) * (cos(heading + kappa(i) * L(i)) - cos(heading));
                            obj.precomp_jac.Xn(i,3+n+j) = obj.precomp_jac.Xn(i-1,3+n+j) + kappa(j)/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));
                            
                            obj.precomp_jac.Yn(i,3+j) = obj.precomp_jac.Yn(i-1,3+j) + 1/kappa(i) * L(j) * (sin(heading + kappa(i) * L(i)) - sin(heading));
                            obj.precomp_jac.Yn(i,3+n+j) = obj.precomp_jac.Yn(i-1,3+n+j) + kappa(j)/kappa(i) * (sin(heading + kappa(i) * L(i)) - sin(heading));
                        else
                            obj.precomp_jac.Xn(i,3+j) = L(i)/kappa(i) * cos(heading + kappa(i) * L(i)) - 1/kappa(i)^2 * (sin(heading + kappa(i) * L(i)) - sin(heading));
                            obj.precomp_jac.Xn(i,3+n+j) = cos(heading + kappa(i) * L(i));

                            obj.precomp_jac.Yn(i,3+j) = L(i)/kappa(i) * sin(heading + kappa(i) * L(i)) + 1/kappa(i)^2 * (cos(heading + kappa(i) * L(i)) - cos(heading));
                            obj.precomp_jac.Yn(i,3+n+j) = sin(heading + kappa(i) * L(i));
                        end
                    end
                end
                
                heading = heading + kappa(i) * L(i);
            end

        end
        
    end

    methods(Static)
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
                
        %% Apply ZOH for association matrix
        function arr = ZOH(arr_in)
            % For all zeros inside the association matrix, change to
            % appropriate close indices
            arr = arr_in;
            met_idx = 0;
            for i=1:length(arr_in)
                if arr_in(i) == 0
                    if met_idx == 0
                        arr(i) = 1;
                    else
                        arr(i) = met_idx;
                    end
                else
                    met_idx = arr_in(i);
                end
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
        
    end
end