classdef ArcFit < handle
    %ArcFit NLS based simple Arc Fitting given data points
    % Given data points and initial arc parameters, this solver module
    % finds the optimal arc parameter set.
    % Since the goal is to find a good 

    properties(Access = private)
        params
        points
        valid = false
        assoc
        options = struct()
    end

    methods(Access = public)        
        %% Constructor
        function obj = ArcFit(params,points)
            obj.params = params;
            obj.points = points;
            obj.assoc = zeros(1,size(obj.points,2));
            
            obj.options.CostThres = 1e-6;
            obj.options.StepThres = 1e-6;
            obj.options.IterThres = 50;
            obj.options.TR = struct();            
            obj.options.TR.eta1 = 0.6;
            obj.options.TR.eta2 = 0.9;
            obj.options.TR.gamma1 = 0.1;
            obj.options.TR.gamma2 = 2;
            obj.options.TR.thres = 1e-6; % Trust Region Radius Threshold

            obj.optimize();
        end
            
        %% Retrieve params 
        function params = getParams(obj)
            params = obj.params;
        end

    end
    methods(Access = private)
        %% Optimize
        function obj = optimize(obj)
            
            while ~obj.valid
                % Perform Data Association
                obj.associate();
                % Trust Region Optimization (Full Convergence)
                obj.TrustRegion();
                % Validate Current Optimization Results
                % If invalid, add more segments
                obj.validate();
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
                        
                        obj.retract(x0,'final');
    
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
                    if tr_rad < obj.options.TRthres
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
            res = vertcat(ME_res,AM_res);
            jac = vertcat(ME_jac,AM_jac);
        end

        %% Create Measurement Block
        function [res,jac] = CreateMEBlock(obj)
            
        end

        %% Create Anchor Measurement Block
        function [res,jac] = CreateAMBlock(obj)
            
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
        end
        
        %% Retract delta values
        function obj = retract(obj,dX)
            n = (length(dX)-3)/2;
            obj.params.x0 = obj.params.x0 + dX(1);
            obj.params.y0 = obj.params.y0 + dX(2);
            obj.params.tau0 = obj.params.tau0 + dX(3);
            obj.params.kappa = obj.params.kappa + dX(4:3+n)';
            obj.params.L = obj.params.L + dX(4+n:end)';
            
            % Perform data association
            obj.associate(); 
        end

        %% Check Optimization Validity
        function obj = validate(obj)
            
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

        %% Find match for data association
        function idx = findMatch(nodePoints,point)
            n = size(nodePoints,2);
            idx = 0;
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
                    idx = i;                    
                end
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

    end
end