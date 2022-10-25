classdef ArcFit < handle
    %ArcFit NLS based simple Arc Fitting given data points
    % Given data points and initial arc parameters, this solver module
    % finds the optimal arc parameter set.
    % Since the goal is to find a good 

    properties(Access = public)
        id
        params
        points
        valid = false
        assoc
        validity
        options = struct()
        opt = struct()
    end

    methods(Access = public)        
        %% Constructor
        function obj = ArcFit(params,points,id)
            obj.id = id;
            obj.params = params;
            obj.points = points(1:2,:);
            obj.assoc = zeros(1,size(obj.points,2));
            
            obj.options.CostThres = 1;
            obj.options.StepThres = 1e-6;
            obj.options.IterThres = inf;
            obj.options.TR = struct();            
            obj.options.TR.eta1 = 0.6;
            obj.options.TR.eta2 = 0.9;
            obj.options.TR.gamma1 = 0.1;
            obj.options.TR.gamma2 = 2;
            obj.options.TR.thres = 1e-6; % Trust Region Radius Threshold
            obj.options.minL = 5; % Minimum arc length constraint
        end
            
        %% Retrieve params 
        function params = getParams(obj)
            params = obj.params;
        end
        
        %% Optimize
        function obj = optimize(obj)
            disp(['Performing Optimization for Segment ID: ',num2str(obj.id)])            
            obj.TrustRegion();
        end
        
        %% Visualize
        function visualize(obj)
            figure(1);
            p_est = plot(obj.points(1,:),obj.points(2,:),'r.'); hold on; grid on; axis equal;
%             plot(obj.points(1,72:78),obj.points(2,72:78),'mx');
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
%             p_seg = plot(SegPoints(1,1),SegPoints(2,1),'ms');
%             plot(SegPoints(1,end),SegPoints(2,end),'ms');
            
            xlabel('Global X(m)'); ylabel('Global Y(m)');
            title('Optimized Vehicle Trajectory and Lane Segments');
            legend([p_est,p_lane,p_node], ...
                   'Input Lane Points', ...
                   'Optimized Arc Spline', ...
                   'Lane Sub-segment');

        end
        
        
    end
    methods(Access = private)  
        %% Create initial parameters
        function obj = initParams(obj)
            P1 = obj.points(:,1); P2 = obj.points(:,end);
            obj.params.tau0 = (P2(2) - P1(2))/(P2(1) - P1(1));
            obj.params.x0 = P1(1); obj.params.y0 = P1(2);
            obj.params.kappa = 1e-6; obj.params.L = sqrt((P1 - P2)'*(P1-P2));
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
                flags = [];
                flags = [flags abs(prev_cost - cost) > obj.options.CostThres];
                flags = [flags step_size > obj.options.StepThres];
                flags = [flags i < obj.options.IterThres];
                
                % Check for oscillation around the local minima
                if length(cost_stack) >= 5
                    osc_flag = obj.DetectOsc(cost_stack);
                else
                    osc_flag = false;
                end

                if length(find(flags)) ~= length(flags) % If any of the criterion is not met, end loop
                    
                    obj.retract(x0);

                    disp('[Optimization Finished...]')
                    idx = find(~flags,1);
                    
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
        function [res,jac] = cost_func(obj,x0)
            obj.retract(x0);
            [ME_res,ME_jac] = obj.CreateMEBlock();
            [AM_res,AM_jac] = obj.CreateAMBlock();
%             [AL_res,AL_jac] = obj.CreateALBlock();
            AL_res = []; AL_jac = [];
            res = vertcat(ME_res,AM_res,AL_res);
            jac = sparse(vertcat(ME_jac,AM_jac,AL_jac));
        end

        %% Create Measurement Block
        function [res,jac] = CreateMEBlock(obj)

            blk_height = nnz(obj.assoc);
            n = length(obj.params.kappa);
            blk_width = 2 * n + 3;
            
            res = zeros(blk_height,1);
            jac = zeros(blk_height,blk_width);
            cnt = 1;
            for i=1:length(obj.points)
                SegIdx = obj.assoc(i);
                Point = obj.points(1:2,i);
                if SegIdx > 0
                    [init_jac,kappa_jac,L_jac,anchored_res] = obj.MEjac(SegIdx,Point);
                    
                    res(cnt) = anchored_res;
                    % Since the number of parameters are not huge, we skip
                    % converting jacobian into sparse array
                    jac(cnt,1:3) = init_jac;
                    jac(cnt,4:3+SegIdx) = kappa_jac;
                    jac(cnt,4+n:3+n+SegIdx) = L_jac;
                    cnt = cnt + 1;
                end
            end

        end
        
        %% Create Measurement Jacobian 
        function [init_jac,kappa_jac,L_jac,anchored_res] = MEjac(obj,SegIdx,Point)
            cov = 0.1^2; eps = 1e-8; eps_kappa = 1e-12;
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa(1:SegIdx); L = obj.params.L(1:SegIdx);
            minL = obj.options.minL;
            anchored_res = obj.MEres(initParams,kappa,L,Point);
            
            init_jac = zeros(1,3);
            kappa_jac = zeros(1,SegIdx);
            L_jac = zeros(1,SegIdx);
            % 1: Perturbation of initial Parameters
            for i=1:3
                initParams_ptb_vec = zeros(1,3);
                initParams_ptb_vec(i) = eps;
                initParams_ptbP = initParams + initParams_ptb_vec;
                initParams_ptbM = initParams - initParams_ptb_vec;
                res_ptbP = obj.MEres(initParams_ptbP,kappa,L,Point);
                res_ptbM = obj.MEres(initParams_ptbM,kappa,L,Point);
                init_jac(i) = (res_ptbP - res_ptbM)/(2*eps);
            end
            init_jac = InvMahalanobis(init_jac,cov);

            % 2: Perturbation of kappa
            for i=1:SegIdx
                kappa_ptb_vec = zeros(1,SegIdx);
                kappa_ptb_vec(i) = eps_kappa;
                kappa_ptbP = kappa + kappa_ptb_vec;
                kappa_ptbM = kappa - kappa_ptb_vec;
                res_ptbP = obj.MEres(initParams,kappa_ptbP,L,Point);
                res_ptbM = obj.MEres(initParams,kappa_ptbM,L,Point);
                kappa_jac(i) = (res_ptbP - res_ptbM)/(2*eps_kappa);
            end
            kappa_jac = InvMahalanobis(kappa_jac,cov);
            
            % 3 Perturbation of L
            for i=1:SegIdx
                L_ptb_vec = zeros(1,SegIdx);
                L_ptb_vec(i) = eps;
                L_ptbP = ((L - minL).^(1/2) + L_ptb_vec).^2 + minL;
                L_ptbM = ((L - minL).^(1/2) - L_ptb_vec).^2 + minL;
                res_ptbP = obj.MEres(initParams,kappa,L_ptbP,Point);
                res_ptbM = obj.MEres(initParams,kappa,L_ptbM,Point);
                L_jac(i) = (res_ptbP - res_ptbM)/(2*eps);
            end
            L_jac = InvMahalanobis(L_jac,cov);

            anchored_res = InvMahalanobis(anchored_res,cov);
        end
        
        %% Create Measurement Residual
        function res = MEres(obj,initParams,kappa,L,Point)
            centerPoints = obj.propCenter(initParams,kappa,L);
            Xc = centerPoints(:,end);

            res = sqrt((Xc - Point)'*(Xc - Point)) - abs(1/kappa(end));
        end

        %% Create Anchor Measurement Block
        function [res,jac] = CreateAMBlock(obj)
            blk_height = 4;
            n = length(obj.params.kappa);
            blk_width = 2 * n + 3;
            res = zeros(blk_height,1);
            jac = zeros(blk_height,blk_width);

            % First
            [init_jac,kappa_jac,L_jac,anchored_res] = obj.AMjac('first');            
            res(1:2) = anchored_res;
            jac(1:2,1:3) = init_jac; 
            jac(1:2,3+1:3+n) = kappa_jac;
            jac(1:2,3+n+1:end) = L_jac;

            % Last
            [init_jac,kappa_jac,L_jac,anchored_res] = obj.AMjac('last');
            res(3:4) = anchored_res;
            jac(3:4,1:3) = init_jac; 
            jac(3:4,3+1:3+n) = kappa_jac;
            jac(3:4,3+n+1:end) = L_jac;
        end

        %% Create Anchor Measurement Jacobian
        function [init_jac,kappa_jac,L_jac,anchored_res] = AMjac(obj,flag)
            cov = diag([1e-5, 1e-5]); eps = 1e-8; eps_kappa = 1e-12;
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            if strcmp(flag,'first')
                kappa = []; L = [];
                anchored_res = obj.AMres(initParams,kappa,L,flag);
                
                kappa_jac = zeros(2,length(obj.params.kappa));
                L_jac = zeros(2,length(obj.params.L));

                init_jac = zeros(2,3);

                for i=1:3
                    initParams_ptb_vec = zeros(1,3);
                    initParams_ptb_vec(i) = eps;
                    initParams_ptbP = initParams + initParams_ptb_vec;
                    initParams_ptbM = initParams - initParams_ptb_vec;
                    res_ptbP = obj.AMres(initParams_ptbP,kappa,L,flag);
                    res_ptbM = obj.AMres(initParams_ptbM,kappa,L,flag);
                    init_jac(:,i) = (res_ptbP - res_ptbM)/(2*eps);
                end
                init_jac = InvMahalanobis(init_jac,cov);
                anchored_res = InvMahalanobis(anchored_res,cov);

            elseif strcmp(flag,'last')
                kappa = obj.params.kappa; L = obj.params.L;
                minL = obj.options.minL;
                anchored_res = obj.AMres(initParams,kappa,L,flag);
                
                kappa_jac = zeros(2,length(obj.params.kappa));
                L_jac = zeros(2,length(obj.params.L));

                init_jac = zeros(2,3);                
                % 1: Perturbation of initial parameters 
                for i=1:3
                    initParams_ptb_vec = zeros(1,3);
                    initParams_ptb_vec(i) = eps;
                    initParams_ptbP = initParams + initParams_ptb_vec;
                    initParams_ptbM = initParams - initParams_ptb_vec;
                    res_ptbP = obj.AMres(initParams_ptbP,kappa,L,flag);
                    res_ptbM = obj.AMres(initParams_ptbM,kappa,L,flag);
                    init_jac(:,i) = (res_ptbP - res_ptbM)/(2*eps);
                end
                init_jac = InvMahalanobis(init_jac,cov);

                for i=1:length(kappa)
                    kappa_ptb_vec = zeros(1,length(kappa));
                    kappa_ptb_vec(i) = eps_kappa;
                    kappa_ptbP = kappa + kappa_ptb_vec;
                    kappa_ptbM = kappa - kappa_ptb_vec;
                    res_ptbP = obj.AMres(initParams,kappa_ptbP,L,flag);
                    res_ptbM = obj.AMres(initParams,kappa_ptbM,L,flag);
                    kappa_jac(:,i) = (res_ptbP - res_ptbM)/(2*eps_kappa);
                end
                kappa_jac = InvMahalanobis(kappa_jac,cov);
                
                % To prevent singularity, after replication, use different
                % epsilon values for jacobian computation
                eps2 = 1e-6 * linspace(1,length(L),length(L));
                for i=1:length(L)
                    L_ptb_vec = zeros(1,length(L));
                    L_ptb_vec(i) = eps2(i);
                    L_ptbP = ((L - minL).^(1/2) + L_ptb_vec).^2 + minL;
                    L_ptbM = ((L - minL).^(1/2) - L_ptb_vec).^2 + minL;
                    res_ptbP = obj.AMres(initParams,kappa,L_ptbP,flag);
                    res_ptbM = obj.AMres(initParams,kappa,L_ptbM,flag);
                    L_jac(:,i) = (res_ptbP - res_ptbM)/(2*eps2(i));
                end
                L_jac = InvMahalanobis(L_jac,cov);

                anchored_res = InvMahalanobis(anchored_res,cov);
            end
        end

        %% Create Anchor Measurement Residual
        function res = AMres(obj,initParams,kappa,L,flag)
            if strcmp(flag,'first')
                X = initParams(1:2)';
                res = X - obj.points(1:2,1);
            elseif strcmp(flag,'last')
                nodePoints = obj.propNode(initParams,kappa,L);
                res = nodePoints(:,end) - obj.points(1:2,end);
            end
        end
        
        %% Create Arc Length Measurement Block
        function [res,jac] = CreateALBlock(obj)
            blk_height = length(obj.params.L);
            n = length(obj.params.kappa);
            blk_width = 2 * n + 3;
            res = zeros(blk_height,1);
            jac = zeros(blk_height,blk_width);

            for i=1:length(obj.params.L)
                res(i) = obj.params.L(i)^(1/2);
                jac(i,3+n+i) = 1/2 * (obj.params.L(i))^(-1/2);
                % 1/(2 * sqrt(obj.params.L(i)));
            end
        end

        %% Data Association
        function obj = associate(obj)
            %ASSOCIATE Data Association 
            % Matches which point belongs to which segment

            % Propagate node points
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa; L = obj.params.L;
            nodePoints = obj.propNode(initParams,kappa,L);
            
            % Find Match
            for i=1:length(obj.points)
                % If not matched, 0 is returned
                obj.assoc(i) = obj.findMatch(nodePoints,obj.points(:,i));
            end

            % Apply zero-order hold for zeros
            obj.assoc = obj.ZOH(obj.assoc);
        end

        %% Retract delta values
        function obj = retract(obj,dX)
            n = (length(dX)-3)/2;
            minL = obj.options.minL;
            obj.params.x0 = obj.params.x0 + dX(1);
            obj.params.y0 = obj.params.y0 + dX(2);
            obj.params.tau0 = obj.params.tau0 + dX(3);
            obj.params.kappa = obj.params.kappa + dX(4:3+n)';
%             obj.params.L = ((obj.params.L).^(1/2) + dX(4+n:end)').^2;
            obj.params.L = ((obj.params.L - minL).^(1/2) + dX(4+n:end)').^2 + minL;
            
%             % Perform data association
            obj.associate(); 
        end

        %% Check Optimization Validity
        function obj = validate(obj)
            obj.validity = zeros(1,length(obj.params.kappa));
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa; L = obj.params.L;
            centerPoints = obj.propCenter(initParams,kappa,L);
            
            % Check normalized Mahalanobis distance error for every data
            % points, and count the number of invalid measurements for each
            % segments
            for i=1:length(obj.points)
                SegIdx = obj.assoc(i);
                if SegIdx > 0
                    if ~obj.IsValid(obj.points(:,i),centerPoints(:,SegIdx),kappa(SegIdx))
                        obj.validity(SegIdx) = obj.validity(SegIdx) + 1;
                    end
                end
            end
            
            if ~isempty(find(obj.validity > 2,1))
                obj.valid = false;
                
                % If current optimization contains invalid segments, find
                % the "worst" segment and divide it into 2.
                
                [~,SegIdx] = max(obj.validity);
                obj.replicate(SegIdx);
                disp(['[Adding Segment at idx ',num2str(SegIdx),']'])
            else
                obj.valid = true;
                disp(['Segment ID: ',num2str(obj.id),' Optimization Finished'])
            end
%             obj.associate();
        end

        %% Replicate current invalid segment
        function obj = replicate(obj,SegIdx)
            if length(obj.params.kappa) == 1
                curr_kappa = obj.params.kappa(SegIdx);
                curr_L = obj.params.L(SegIdx);
                obj.params.kappa = [curr_kappa, curr_kappa];
                obj.params.L = [1/2 * curr_L, 1/2 * curr_L];
            else
                if SegIdx == 1
                    curr_kappa = obj.params.kappa(SegIdx);    
                    next_kappa = obj.params.kappa(SegIdx+1);
    %                 new_kappa = 2/(1/curr_kappa + 1/next_kappa);
                    new_kappa = curr_kappa;
%                     new_kappa = 1/2 * (curr_kappa + next_kappa);
                    curr_L = obj.params.L(SegIdx);
                    next_L = obj.params.L(SegIdx+1);
                    new_L = 1/3 * (curr_L + next_L);
                    obj.params.kappa = [curr_kappa, new_kappa, obj.params.kappa(2:end)];
                    
%                     obj.params.L = [1/2 * curr_L, 1/2 * curr_L, obj.params.L(2:end)];
                    obj.params.L = [new_L, new_L, new_L obj.params.L(3:end)];
                elseif SegIdx == length(obj.params.kappa)
                    curr_kappa = obj.params.kappa(SegIdx);
%                     prev_kappa = obj.params.kappa(SegIdx-1);
%                     new_kappa = 1/2 * (prev_kappa + curr_kappa);
    %                 new_kappa = 2/(1/prev_kappa + 1/curr_kappa);
                    new_kappa = curr_kappa; 
                    curr_L = obj.params.L(SegIdx);
                    prev_L = obj.params.L(SegIdx-1);
                    new_L = 1/3 * (curr_L + prev_L);
    
                    obj.params.kappa = [obj.params.kappa(1:end-1), new_kappa, curr_kappa];
                    obj.params.L = [obj.params.L(1:end-2), new_L, new_L, new_L];
%                     obj.params.L = [obj.params.L(1:end-1), 1/2 * curr_L, 1/2 * curr_L];
                else
                    curr_kappa = obj.params.kappa(SegIdx);
                    front_kappa = obj.params.kappa(1:SegIdx-1);
                    back_kappa = obj.params.kappa(SegIdx+1:end);
                    
    %                 front_new_kappa = 2/(1/curr_kappa + 1/front_kappa(end));
    %                 back_new_kappa = 2/(1/curr_kappa + 1/back_kappa(1));
                    front_new_kappa = curr_kappa;
                    back_new_kappa = curr_kappa;
%                     front_new_kappa = 1/2 * (curr_kappa + front_kappa(end));
%                     back_new_kappa = 1/2 * (curr_kappa + back_kappa(1));
    
                    curr_L = obj.params.L(SegIdx);
                    front_L = obj.params.L(1:SegIdx-1);
                    back_L = obj.params.L(SegIdx+1:end);
                    
                    new_L = 1/4 * (curr_L + front_L(end) + back_L(1));                
    
                    obj.params.kappa = [front_kappa, front_new_kappa, back_new_kappa, back_kappa];
                    obj.params.L = [front_L(1:end-1), new_L, new_L, new_L, new_L, back_L(2:end)];
%                     obj.params.L = [front_L, 1/2 * curr_L, 1/2 * curr_L, back_L];
                end
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
%             disp(centerPoints)
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

        %% Find match for data association
        function idx = findMatch(nodePoints,point)
            n = size(nodePoints,2);
            cand_idx = [];
            VertDist = [];
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
                    VertDist = [VertDist sqrt((X_vert - point)' * (X_vert - point))];
                end
            end

            % If there are many possible matches, compare each match
            % Using the vertical distance
            
            if isempty(cand_idx)
                idx = 0;
            else
                [~,min_idx] = min(VertDist);
                idx = cand_idx(min_idx);
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
        
        %% Return validity of certain point after optimization
        function flag = IsValid(X,Xc,kappa)
            cov = 0.5^2;
            diff = sqrt((X - Xc)' * (X - Xc)) - abs(1/kappa);
            Ndist = diff' / cov * diff;
            if Ndist > chi2inv(0.99,2)
                flag = false;
            else
                flag = true;
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

    end
end