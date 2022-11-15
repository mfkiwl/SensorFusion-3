classdef ArcFit2 < handle
    properties(Access = public)
        id
        params        
        points  
        est_points = {}
        covs        
        assoc
        fullL
        blk_width 
        ref_idxs
        opt
        nodePoints
    end

    methods (Access = public)
        %% Constructor
        function obj = ArcFit2(params,points,covs,id)
            obj.id = id;
            obj.params = params;            
            obj.points = points;            
            obj.covs = covs;
            
            % Association is fixed throughout the optimization process
            obj.associate();
%             obj.associateR();

            % Optimization conditions
            obj.opt.options = struct();
            obj.opt.options.CostThres = 1e-6;
            obj.opt.options.StepThres = 1e-6;
            obj.opt.options.IterThres = 1.5e3;
            obj.opt.options.Algorithm = 'TR';
            obj.opt.options.TR = struct();
            obj.opt.options.TR.thres = 1e-6;
            obj.opt.options.TR.eta1 = 0.5;
            obj.opt.options.TR.eta2 = 0.9;
            obj.opt.options.TR.gamma1 = 0.1;  
            obj.opt.options.TR.gamma2 = 2;
        end
        
        %% Retrieve params 
        function params = getParams(obj)
            params = obj.params;
        end
        
        %% Optimizer
        function obj = optimize(obj)
            X0 = obj.createInitVal();
            obj.blk_width = length(X0);
            obj.opt.x0 = zeros(obj.blk_width,1);
            obj.TrustRegion();
%             jac_p = sparse([0,0,0,0,1,0]);
%             options = optimoptions('lsqnonlin', ...                                                                      
%                                    'MaxFunctionEvaluations',5e3, ...
%                                    'MaxIterations',5e3, ...
%                                    'Display','iter-detailed', ...
%                                    'SpecifyObjectiveGradient',true, ...
%                                    'CheckGradients',false);
%             disp(['[Performing Optimization for Segment ID: ',num2str(obj.id),']']) 
% 
%             % Need to add validation and replication functions
%             n = length(obj.params.kappa);
%             lb = [repmat(-inf,1,3+n),zeros(1,size(obj.points,2)-1)];
%             X = lsqnonlin(@obj.cost_func,X0,[],[],options);
            
            % arc length variable should be strictly increasing : check

            
            % Save Optimization Results
%             n = length(obj.params.kappa);
%             obj.params.x0 = X(1);
%             obj.params.y0 = X(2);
%             obj.params.tau0 = X(3);
%             obj.params.kappa = X(3+1:3+n)';
%             
%             obj.fullL = [0,X(3+n+1:end)'];
%             obj.params.L = obj.fullL(obj.ref_idxs);
            
            % Propagate Points and save optimziation results
            obj.propPoints();
        end

        %% Visualize
        function visualize(obj)
            figure(30); hold on; grid on; axis equal;
            plot(obj.points(1,:),obj.points(2,:),'r.');

            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa;
            L = obj.params.L;
            for i=1:length(obj.est_points)
                plot(obj.est_points{i}(1,:),obj.est_points{i}(2,:),'b.');
                filledL = linspace(0,obj.fullL(obj.ref_idxs(i)),1e3);
                PC = [];
                
                for j=1:length(filledL)
                    PC = [PC, obj.getPoint(i,initParams,kappa,L,filledL(j))];
                end
                plot([PC(1,1),PC(1,end)],[PC(2,1),PC(2,end)],'ms');
                plot(PC(1,:),PC(2,:),'k--');
            end
            
            
        end

    end

    methods (Access = private)
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
                
%                 if ~isempty(find(isnan(h_gn),1)) || ~isempty(find(isnan(h_gd),1))
%                     obj.opt.info = A;
%                     obj.errorAnalysis();
%                 end

                x0 = obj.ComputeDogLeg(h_gn,h_gd,tr_rad);
                
                dummy_params = obj.params;
                dummy_fullL = obj.fullL;
                
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
                    obj.params = dummy_params;                                   
                    obj.fullL = dummy_fullL;
                    
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
                        
                        obj.retract(x0);
    
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
                        obj.retract(x0);
                        disp(['Current trust region radius ',num2str(tr_rad),' is below threshold: ',num2str(obj.opt.options.TR.thres)])
                        break;
                    end
                end
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
                
                osc_flag = false;
%                 % Check for oscillation around the local minima
%                 if length(cost_stack) >= 5
%                     osc_flag = obj.DetectOsc(cost_stack);
%                 else
%                     osc_flag = false;
%                 end

                if length(find(obj.opt.flags)) ~= length(obj.opt.flags) % If any of the criterion is not met, end loop
                    
                    obj.retract(x0);

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

        %% NLS Cost Function
        function [res,jac] = cost_func(obj,X0)
            
            obj.retract(X0);

%             n = length(obj.params.kappa);
%             initParams = X0(1:3);
%             kappa = X0(3+1:3+n)';
%             Lf = [0, X0(3+n+1:end)']; % Arc length for all data points
%             % Be aware that the first 0 in Lf is not to be optimized
%             % It is added just for computation simplicity
%             
%             L = Lf(obj.ref_idxs); % Node Point Lengths            


            % Pre compute Jacobians
%             [precompX, precompY] = obj.precomputeJac(initParams,kappa,L);
            
            [AL_res,AL_jac] = obj.createALBlock();
            [ME_res,ME_jac] = obj.createMEBlock();
            [AC_res,AC_jac] = obj.createACBlock();

            res = vertcat(AL_res,ME_res,AC_res);
            jac = vertcat(AL_jac,ME_jac,AC_jac);
            obj.opt.jac = jac;
        end
        
        %% Retract Variables
        function obj = retract(obj,delta)
            n = length(obj.params.kappa);
            x0_delta = delta(1);
            y0_delta = delta(2);
            tau0_delta = delta(3);
            kappa_delta = delta(3+1:3+n)';
            sqrtL_delta = delta(3+n+1:end)';
            
            obj.params.x0 = obj.params.x0 + x0_delta;
            obj.params.y0 = obj.params.y0 + y0_delta;
            obj.params.tau = obj.params.tau0 + tau0_delta;
            obj.params.kappa = obj.params.kappa + kappa_delta;
            obj.fullL = (obj.fullL.^(1/2) + [0, sqrtL_delta]).^2;
            obj.params.L = obj.fullL(obj.ref_idxs);

%             obj.nodePoints = obj.propNode([obj.params.x0,obj.params.y0,obj.params.tau0], ...
%                                            obj.params.kappa,obj.params.L);
        end
        
        %% Create Arc Length Block
        function [res,jac] = createALBlock(obj)
            blk_height = size(obj.points,2)-1;
            res = zeros(blk_height,1);
            bnds = obj.params.bnds;
            n = length(obj.params.kappa); % careful for replication
            cnt = 1;
            jac_cnt = 0;
            I = zeros(1,2*blk_height-size(bnds,1)); J = I; V = I;
            for i=1:size(bnds,1)
                lb = bnds(i,1); ub = bnds(i,2);
                
                prevL = 0;
                for j=lb+1:ub                    
                    currL = obj.fullL(j);
%                     if currL < prevL
%                         disp([i,j])
%                         error('Arc Length should be strictly increasing')
%                     end

                    dist = norm(obj.points(:,j) - obj.points(:,j-1));
                    if dist > 10 
                        % If points are far away from each other, straight
                        % line distance will be quite different from arc
                        % length distance.
                        cov = 1^2;
                    else
                        cov = 0.1^2;
                    end
                    res(cnt) = InvMahalanobis(currL - prevL - dist,cov);
                    
                    if j ~= lb+1
                        [I1,J1,V1] = sparseFormat(cnt,3+n+j-2:3+n+j-1,InvMahalanobis([-2*sqrt(prevL),2*sqrt(currL)],cov));
                        I(jac_cnt+1:jac_cnt+2) = I1;
                        J(jac_cnt+1:jac_cnt+2) = J1;
                        V(jac_cnt+1:jac_cnt+2) = V1;
                        jac_cnt = jac_cnt + 2;
                    else
                        [I1,J1,V1] = sparseFormat(cnt,3+n+j-1,InvMahalanobis(2*sqrt(currL),cov));
                        I(jac_cnt+1) = I1;
                        J(jac_cnt+1) = J1;
                        V(jac_cnt+1) = V1;
                        jac_cnt = jac_cnt + 1;
                    end

                    cnt = cnt + 1;
                    prevL = currL;
                end
            end
            jac = sparse(I,J,V,blk_height,obj.blk_width);

        end

        %% Create Measurement Block
        function [res,jac] = createMEBlock(obj)
            n = length(obj.params.kappa);
            m = size(obj.points,2)-1;
%             m = size(obj.points,2);
            blk_height = 2*m;
%             blk_heightR = 2*size(obj.points,2);

%             heading = initParams(3);
            
            res = zeros(blk_height,1);
%             resR = zeros(blk_heightR,1);
            cnt = 1; 
%             cntR = 1;
%             jac_cnt = 0; 
%             jac_cntR = 0;

%             I = zeros(1,2*(2*sum(obj.assoc(2:end))+3*m)); J = I; V = I;
%             IR = zeros(1,2*(2*sum(obj.assoc)+3*m)); JR = IR; VR = IR;
            I = []; J = []; V = [];
            for i=1:n
                Idxs = find(obj.assoc == i);
%                 IdxsR = find(obj.assoc == i);
                
                               

                % 0m Points
                if i == 1
                    jstart = 2;
                else
                    jstart = 1;
                end
                for j=jstart:length(Idxs)
                    
%                     cov = diag([1e-2,1e-2]);
                    cov = reshape(obj.covs(:,Idxs(j)),2,2);
                    if i == n && j == length(Idxs)
                        cov = diag([1e-4,1e-4]);
                    end
                    [init_jac,kappa_jac,L_jac,anchored_res] = obj.getMEJac(Idxs(j),i);

%                     headingN = heading + kappa(i) * L_;
%                     
%                     % Residual Computation
%                     xij = xim1 + 1/kappa(i) * (sin(headingN) - sin(heading));
%                     yij = yim1 - 1/kappa(i) * (cos(headingN) - cos(heading));
% 
                    res(2*cnt-1:2*cnt) = InvMahalanobis(anchored_res,cov);
%                     
%                     % Measurement Jacobian Computation
%                     init_jac = zeros(2,3);
%                     init_jac(1,1) = 1; init_jac(2,2) = 1;
% 
%                     kappa_jac = zeros(2,i);
%                     L_jac = zeros(2,i);
                    
%                     if i == 1 % Matched with first segment
%                         init_jac(1,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
%                         init_jac(2,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
%                         
%                         % k = 1 only
%                         for k=1:i
%                             kappa_jac(1,k) = L_/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                             kappa_jac(2,k) = L_/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                             
%                             L_jac(1,k) = cos(headingN);
%                             L_jac(2,k) = sin(headingN);
%                         end
%                     else
%                         init_jac(1,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
%                         init_jac(2,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));
%                         
%                         for k=1:i
%                             if k ~= i
%                                 kappa_jac(1,k) = precompX(i-1,3+k) + L(k)/kappa(i) * (cos(headingN) - cos(heading));
%                                 kappa_jac(2,k) = precompY(i-1,3+k) + L(k)/kappa(i) * (sin(headingN) - sin(heading));
% 
%                                 L_jac(1,k) = precompX(i-1,3+n+k) + kappa(k)/kappa(i) * (cos(headingN) - cos(heading));
%                                 L_jac(2,k) = precompY(i-1,3+n+k) + kappa(k)/kappa(i) * (sin(headingN) - sin(heading));
%                             else
%                                 kappa_jac(1,k) = L_/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                                 kappa_jac(2,k) = L_/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                                 
%                                 L_jac(1,k) = cos(headingN);
%                                 L_jac(2,k) = sin(headingN);
%                             end
%                         end
%                     end

                    % Normalize
                    init_jac = InvMahalanobis(init_jac,cov);
                    kappa_jac = InvMahalanobis(kappa_jac,cov);
                    L_jac = InvMahalanobis(L_jac,cov);
% 
                    % Assign
                    [I_i,J_i,V_i] = sparseFormat(2*cnt-1:2*cnt,1:3,init_jac);
                    [I_k,J_k,V_k] = sparseFormat(2*cnt-1:2*cnt,3+1:3+i,kappa_jac);
                    [I_sL,J_sL,V_sL] = sparseFormat(2*cnt-1:2*cnt,3+n+[obj.ref_idxs(1:i-1),Idxs(j)]-1,L_jac);
                    
%                     I(jac_cnt+1:jac_cnt+2*(2*i+3)) = [I_i,I_k,I_sL];
%                     J(jac_cnt+1:jac_cnt+2*(2*i+3)) = [J_i,J_k,J_sL];
%                     V(jac_cnt+1:jac_cnt+2*(2*i+3)) = [V_i,V_k,V_sL];
                    I = [I, I_i, I_k, I_sL];
                    J = [J, J_i, J_k, J_sL];
                    V = [V, V_i, V_k, V_sL];

                    cnt = cnt + 1;
%                     jac_cnt = jac_cnt + 2*(2*i+3);
                end
                jac = sparse(I,J,V,blk_height,obj.blk_width);

                % Rest
%                 for j=1:length(IdxsR)
%                     sqrtL_ = sqrtL(IdxsR(j));
%                     cov = reshape(obj.covs(:,IdxsR(j)),2,2);
% 
%                     headingN = heading + kappa(i) * sqrtL_^2;
%                     
%                     % Residual Computation
%                     xij = xim1 + 1/kappa(i) * (sin(headingN) - sin(heading));
%                     yij = yim1 - 1/kappa(i) * (cos(headingN) - cos(heading));
% 
%                     resR(2*cntR-1:2*cntR) = InvMahalanobis([xij;yij] - obj.points(:,IdxsR(j)),cov);
                    
%                     % Measurement Jacobian Computation
%                     init_jacR = zeros(2,3);
%                     init_jacR(1,1) = 1; init_jacR(2,2) = 1;
% 
%                     kappa_jacR = zeros(2,i);
%                     sqrtL_jacR = zeros(2,i);
%                     
%                     if i == 1 % Matched with first segment
%                         init_jacR(1,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
%                         init_jacR(2,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
%                         
%                         % k = 1 only
%                         for k=1:i
%                             kappa_jacR(1,k) = sqrtL_^2/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                             kappa_jacR(2,k) = sqrtL_^2/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                             
%                             sqrtL_jacR(1,k) = cos(headingN) * 2 * sqrtL_;
%                             sqrtL_jacR(2,k) = sin(headingN) * 2 * sqrtL_;
%                         end
%                     else
%                         init_jacR(1,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
%                         init_jacR(2,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));
%                         
%                         for k=1:i
%                             if k ~= i
%                                 kappa_jacR(1,k) = precompX(i-1,3+k) + L(k)/kappa(i) * (cos(headingN) - cos(heading));
%                                 kappa_jacR(2,k) = precompY(i-1,3+k) + L(k)/kappa(i) * (sin(headingN) - sin(heading));
% 
%                                 sqrtL_jacR(1,k) = precompX(i-1,3+n+k) + kappa(k)/kappa(i) * (cos(headingN) - cos(heading)) * 2 * sqrtL_;
%                                 sqrtL_jacR(2,k) = precompY(i-1,3+n+k) + kappa(k)/kappa(i) * (sin(headingN) - sin(heading)) * 2 * sqrtL_;
%                             else
%                                 kappa_jacR(1,k) = sqrtL_^2/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                                 kappa_jacR(2,k) = sqrtL_^2/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                                 
%                                 sqrtL_jacR(1,k) = cos(headingN) * 2 * sqrtL_;
%                                 sqrtL_jacR(2,k) = sin(headingN) * 2 * sqrtL_;
%                             end
%                         end
%                     end
% 
%                     % Normalize
%                     init_jacR = InvMahalanobis(init_jacR,cov);
%                     kappa_jacR = InvMahalanobis(kappa_jacR,cov);
%                     sqrtL_jacR = InvMahalanobis(sqrtL_jacR,cov);
% 
%                     % Assign
%                     [IR_i,JR_i,VR_i] = sparseFormat(2*cntR-1:2*cntR,1:3,init_jacR);
%                     [IR_k,JR_k,VR_k] = sparseFormat(2*cntR-1:2*cntR,3+1:3+i,kappa_jacR);
%                     [IR_sL,JR_sL,VR_sL] = sparseFormat(2*cntR-1:2*cntR,3+n+[obj.ref_idxs(1:i-1),m0+IdxsR(j)],sqrtL_jacR);
%                     
%                     IR(jac_cntR+1:jac_cntR+2*(2*i+3)) = [IR_i,IR_k,IR_sL];
%                     JR(jac_cntR+1:jac_cntR+2*(2*i+3)) = [JR_i,JR_k,JR_sL];
%                     VR(jac_cntR+1:jac_cntR+2*(2*i+3)) = [VR_i,VR_k,VR_sL];

%                     cntR = cntR + 1;
%                     jac_cntR = jac_cntR + 2*(2*i+3);
%                 end
% 
%                 heading = heading + kappa(i) * L(i);
            end

%             jac0 = sparse(I0,J0,V0,blk_height0,obj.blk_width);
%             jacR = sparse(IR,JR,VR,blk_heightR,obj.blk_width);
%             res = vertcat(res0,resR);
%             jac = vertcat(jac0,jacR);

        end

        %% Compute Numerical Jacobian for Measurement Block
        function [init_jac,kappa_jac,L_jac,anchored_res] = getMEJac(obj,PointIdx,SegIdx)
            
            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa;
            L = obj.params.L;
            L_ = obj.fullL(PointIdx);

            anchored_res = obj.getMERes(PointIdx,SegIdx,initParams,kappa,L,L_);
            
            epsilon = 1e-6;
            epsilon_kappa = 1e-10;
            % init_jac
            init_jac = zeros(2,3);
            for i=1:3
                initParamsPtb_vec = zeros(1,3);
                initParamsPtb_vec(i) = epsilon;
                initParamsPtb = initParams + initParamsPtb_vec;
                init_jac(:,i) = (obj.getMERes(PointIdx,SegIdx,initParamsPtb,kappa,L,L_) - anchored_res)/epsilon;
            end
            
            % kappa_jac
            kappa_jac = zeros(2,SegIdx);
            for i=1:SegIdx
                kappaPtb_vec = zeros(1,length(kappa));
                kappaPtb_vec(i) = epsilon_kappa;
                kappaPtb = kappa + kappaPtb_vec;
                kappa_jac(:,i) = (obj.getMERes(PointIdx,SegIdx,initParams,kappaPtb,L,L_) - anchored_res)/epsilon_kappa;
            end

            % L_jac
            L_jac = zeros(2,SegIdx);
            for i=1:SegIdx-1
                LPtb_vec = zeros(1,length(L));
                LPtb_vec(i) = epsilon;
                LPtb = (L.^(1/2) + LPtb_vec).^2;
                L_jac(:,i) = (obj.getMERes(PointIdx,SegIdx,initParams,kappa,LPtb,L_) - anchored_res)/epsilon; 
            end
            
            L_Ptb_vec = epsilon;
            L_Ptb = (L_.^(1/2) + L_Ptb_vec).^2;
            L_jac(:,SegIdx) = (obj.getMERes(PointIdx,SegIdx,initParams,kappa,L,L_Ptb) - anchored_res)/epsilon;

        end
        
        %% Compute Measurement Residual
        function res = getMERes(obj,PointIdx,SegIdx,initParams,kappa,L,L_)
            Point = obj.points(:,PointIdx);
            nodePoints_ = obj.propNode(initParams,kappa(1:SegIdx-1),L(1:SegIdx-1));
            heading = initParams(3) + kappa(1:SegIdx-1) * L(1:SegIdx-1)';
            headingN = heading + kappa(SegIdx) * L_;

            predPoint = nodePoints_(:,end) + 1/kappa(SegIdx) * [sin(headingN) - sin(heading);-cos(headingN) + cos(heading)];

            res = predPoint - Point;
        end

        %% Create Init Point Anchoring Block
        function [res,jac] = createACBlock(obj)
            X = obj.points(:,1);
%             cov = reshape(obj.covs0(:,1),2,2);
            cov = diag([1e-4,1e-4]);

            res = InvMahalanobis(X - [obj.params.x0;obj.params.y0],cov);
            jac = sparse([],[],[],2,obj.blk_width);
            jac(1,1) = -1; jac(2,2) = -1;
            jac = InvMahalanobis(jac,cov);
        end

        %% Data Association 0m Previewed Measurements
        function obj = associate(obj)
            obj.assoc = zeros(1,size(obj.points,2));
            % If state idx is used, then association is fixed throughout
            % the whole optimization process
            bnds = obj.params.bnds;
            for i=1:size(bnds,1)
                lb = bnds(i,1); ub = bnds(i,2);
                if i~=1
                    obj.assoc(lb+1:ub) = i;
                else
                    obj.assoc(lb:ub) = i;
                end
            end

            % Create reference idxs for (end) node points
            obj.ref_idxs = obj.params.bnds(:,2)';
        end

        %% Data Association Remaining Measurements
        function obj = associateR(obj)
%             %ASSOCIATE Data Association 
%             % Matches which point belongs to which segment
%             % This association is for remaining lane points
%             % Not to be used in this research.
%             % For Future Use!
%             obj.assoc = zeros(1,size(obj.points,2));            
%             for i=1:size(obj.points,2)
%                 Px = obj.points(1,i); Py = obj.points(2,i);
%                 P_refx = obj.points0(1,:); P_refy = obj.points0(2,:);
%                 D = sqrt((P_refx - Px).^2 + (P_refy - Py).^2);
%                 [~,idx] = min(D);
%                 
%                 if idx ~= 1 && idx ~= size(obj.points0,2)
%                     % Check if idx belongs to any of the boundary index
%                     [r,~] = find(obj.params.bnds == idx);
%                     
%                     % If nearest 0m point is matched, then determine
%                     % association using nearby 0m points
%                     if ~isempty(r)
%                         prev = obj.points0(idx-1);
%                         next = obj.points0(idx+1);
%                         
%                         d1 = norm([Px;Py] - prev);
%                         d2 = norm([Px;Py] - next);
% 
%                         if d1 >= d2 % Matched to "Next"
%                             obj.assoc(i) = obj.assoc0(idx+1);
%                         else % Matched to "Prev"
%                             obj.assoc(i) = obj.assoc0(idx-1);
%                         end
%                     else
%                         obj.assoc(i) = obj.assoc0(idx);
%                     end
%                 else
%                     obj.assoc(i) = obj.assoc0(idx);
%                 end
%             end                
        end
        
        %% Create Initial Value
        function X0 = createInitVal(obj)
            initParams = [obj.params.x0; obj.params.y0; obj.params.tau0];
            kappa = obj.params.kappa';
%             nodeSqrtL = (obj.params.L).^(1/2);
% 
%             % L: 0m L + Remaining L
%             % Use Association information to assign initial value
%             sqrtL0 = zeros(size(obj.points0,2),1);
%             sqrtL = zeros(size(obj.points,2),1);
%             
%             for i=1:size(obj.points0,2)
%                 SegIdx = obj.assoc0(i);
%                 if i == obj.params.bnds(SegIdx,2)
%                     % Initial value for node points
%                     % sqrt(L)
%                     sqrtL0(i) = nodeSqrtL(SegIdx);
%                 else
%                     % Initial value for other points
%                     % sqrt(L/2)
%                     sqrtL0(i) = sqrt(nodeSqrtL(SegIdx)^2/2);
%                 end
%             end
%             
%             for i=1:size(obj.points,2)
%                 SegIdx = obj.assoc(i);
%                 sqrtL(i) = sqrt(nodeSqrtL(SegIdx)^2/2);
%             end
%             
%             X0 = [initParams;kappa;sqrtL0;sqrtL];

            Ls = [];
            for i=1:size(obj.params.bnds,1)                
                % Lower bound arc length is not a variable!
                lb = obj.params.bnds(i,1);
                ub = obj.params.bnds(i,2);
                L = 0;
                prev_point = obj.points(:,lb);
                for j=lb:ub-1
                    curr_point = obj.points(:,j+1);
                    delL = norm(curr_point - prev_point);
                    Ls = [Ls; L + delL];                    
                    L = L + delL;
                    prev_point = curr_point;
                end
            end
%             Ls = 5*ones(size(obj.points,2)-1,1);
            obj.fullL = [0, Ls'];
            X0 = [initParams;kappa;Ls];
        end
        
        %% Propagate Points after optimization
        function obj = propPoints(obj)
            obj.est_points = {};

            initParams = [obj.params.x0, obj.params.y0, obj.params.tau0];
            kappa = obj.params.kappa;
            L = obj.params.L;
            for i=1:size(obj.params.bnds,1)
                lb = obj.params.bnds(i,1);
                ub = obj.params.bnds(i,2);
                PC = [];
                if i == 1
                    lb_ = lb;
                else
                    lb_ = lb+1;
                end
                for j=lb_:ub
                    L_ = obj.fullL(j);
                    PC = [PC obj.getPoint(i,initParams,kappa,L,L_)];
                end
                obj.est_points = [obj.est_points {PC}];
            end
        end

        %% Compute points using parameters
        function PC = getPoint(obj,SegIdx,initParams,kappa,L,L_)
            Xn = obj.propNode(initParams,kappa(1:SegIdx-1),L(1:SegIdx-1));
            heading = initParams(3) + kappa(1:SegIdx-1) * L(1:SegIdx-1)';
            headingN = heading + kappa(SegIdx) * L_;

            PC = Xn(:,end) + 1/kappa(SegIdx) * [sin(headingN) - sin(heading); -cos(headingN) + cos(heading)];
        end

    end

    methods (Static)              
        %% Propagate Node Points
        function nodePoints = propNode(initParams,kappa,L)
            x0 = initParams(1); y0 = initParams(2); heading = initParams(3);            
            nodePoints = [x0; y0];
            
            for i=1:length(kappa)
                nodePoints = [nodePoints nodePoints(:,end) + 1/kappa(i) * [sin(heading + kappa(i) * L(i)) - sin(heading);
                                                                           -cos(heading + kappa(i) * L(i)) + cos(heading)]];
                heading = heading + kappa(i) * L(i);
            end
        end
        
        %% Precompute Node Points Jacobians
        function [precompX, precompY] = precomputeJac(initParams,kappa,L)
            % L is the arc length of each segment
            % Not state variables
            n = length(kappa);            
            heading = initParams(3);

            precompX = zeros(n,2*n+3); precompY = precompX;            
            
            % Optimization L = s^2 >= 0
%             for i=1:n
%                 headingN = heading + kappa(i) * L(i);
%                 if i == 1
%                     precompX(i,1) = 1;
%                     precompX(i,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
%                     precompX(i,3+i) = L(i)/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                     precompX(i,3+n+i) = cos(headingN) * 2 * sqrtL(i);
% 
%                     precompY(i,2) = 1;
%                     precompY(i,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
%                     precompY(i,3+i) = L(i)/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                     precompY(i,3+n+i) = sin(headingN) * 2 * sqrtL(i);
%                 else                    
%                     precompX(i,1) = precompX(i-1,1);
%                     precompX(i,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
%                     
%                     precompY(i,2) = precompY(i-1,2);
%                     precompY(i,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));
% 
%                     for j=1:i
%                         if j ~= i
%                             precompX(i,3+j) = precompX(i-1,3+j) + L(j)/kappa(i) * (cos(headingN) - cos(heading));
%                             precompX(i,3+n+j) = precompX(i-1,3+n+j) + kappa(j)/kappa(i) * (cos(headingN) - cos(heading)) * 2 * sqrtL(j);
%                             
%                             precompY(i,3+j) = precompY(i-1,3+j) + L(j)/kappa(i) * (sin(headingN) - sin(heading));
%                             precompY(i,3+n+j) = precompY(i-1,3+n+j) + kappa(j)/kappa(i) * (sin(headingN) - sin(heading)) * 2 * sqrtL(j);
%                         else
%                             precompX(i,3+i) = L(i)/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                             precompX(i,3+n+i) = cos(headingN) * 2 * sqrtL(i);
% 
%                             precompY(i,3+i) = L(i)/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                             precompY(i,3+n+i) = sin(headingN) * 2 * sqrtL(i);                        
%                         end
%                     end                        
%                 end
%                 
%                 heading = headingN;
%             end
            
            % Optimization simply with L
            for i=1:n
                headingN = heading + kappa(i) * L(i);
                if i == 1
                    precompX(i,1) = 1;
                    precompX(i,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
                    precompX(i,3+i) = L(i)/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
                    precompX(i,3+n+i) = cos(headingN);

                    precompY(i,2) = 1;
                    precompY(i,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
                    precompY(i,3+i) = L(i)/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
                    precompY(i,3+n+i) = sin(headingN);
                else                    
                    precompX(i,1) = precompX(i-1,1);
                    precompX(i,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
                    
                    precompY(i,2) = precompY(i-1,2);
                    precompY(i,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));

                    for j=1:i
                        if j ~= i
                            precompX(i,3+j) = precompX(i-1,3+j) + L(j)/kappa(i) * (cos(headingN) - cos(heading));
                            precompX(i,3+n+j) = precompX(i-1,3+n+j) + kappa(j)/kappa(i) * (cos(headingN) - cos(heading));
                            
                            precompY(i,3+j) = precompY(i-1,3+j) + L(j)/kappa(i) * (sin(headingN) - sin(heading));
                            precompY(i,3+n+j) = precompY(i-1,3+n+j) + kappa(j)/kappa(i) * (sin(headingN) - sin(heading));
                        else
                            precompX(i,3+i) = L(i)/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
                            precompX(i,3+n+i) = cos(headingN);

                            precompY(i,3+i) = L(i)/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
                            precompY(i,3+n+i) = sin(headingN);                        
                        end
                    end                        
                end
                
                heading = headingN;
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

        %% Detect oscillation by computing maximum deviation from average (recent 5 costs)
        function osc_flag = DetectOsc(cost_stack)
            last_five = cost_stack(end-4:end);
            avg = mean(last_five);
            delta = last_five - avg;
            % If all recent costs are near the average, detect oscillation
            if max(delta) < 1e0
                osc_flag = true;
            else
                osc_flag = false;
            end
        end

    end
end