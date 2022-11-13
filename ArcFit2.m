classdef ArcFit2 < handle
    properties(Access = public)
        id
        params        
        points        
        covs        
        assoc
        fullL
        blk_width 
        ref_idxs
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
            obj.associateR();
        end
        
        %% Retrieve params 
        function params = getParams(obj)
            params = obj.params;
        end
        
        %% Optimizer
        function obj = optimize(obj)
            X0 = obj.createInitVal();
            obj.blk_width = length(X0);
%             jac_p = sparse([0,0,0,0,1,0]);
            options = optimoptions('lsqnonlin', ...                                                                      
                                   'MaxFunctionEvaluations',2000, ...
                                   'MaxIterations',2000, ...
                                   'Display','iter-detailed', ...
                                   'SpecifyObjectiveGradient',true);
            disp(['[Performing Optimization for Segment ID: ',num2str(obj.id),']']) 

            % Need to add validation and replication functions
            n = length(obj.params.kappa);
            lb = [repmat(-inf,1,3+n),zeros(1,size(obj.points,2)-1)];
            X = lsqnonlin(@obj.cost_func,X0,lb,[],options);
            
            % arc length variable should be strictly increasing : check

            
            % Save Optimization Results
            n = length(obj.params.kappa);
            obj.params.x0 = X(1);
            obj.params.y0 = X(2);
            obj.params.tau0 = X(3);
            obj.params.kappa = X(3+1:3+n)';
            
            obj.fullL = [0,X(3+n+1:end)'];
            obj.params.L = obj.fullL(obj.ref_idxs);
            
        end

        %% Visualize
        function visualize(obj)

        end

    end

    methods (Access = private)
        %% NLS Cost Function
        function [res,jac] = cost_func(obj,X0)
            n = length(obj.params.kappa);
            initParams = X0(1:3);
            kappa = X0(3+1:3+n)';
            Lf = [0, X0(3+n+1:end)']; % Arc length for all data points
            % Be aware that the first 0 in Lf is not to be optimized
            % It is added just for computation simplicity
            
            L = Lf(obj.ref_idxs); % Node Point Lengths            

            nodePoints = obj.propNode(initParams,kappa,L);
            % Pre compute Jacobians
            [precompX, precompY] = obj.precomputeJac(initParams,kappa,L);
            
            [AL_res,AL_jac] = obj.createALBlock(Lf);
            [ME_res,ME_jac] = obj.createMEBlock(initParams,kappa,Lf,L,nodePoints,precompX,precompY);
            [AC_res,AC_jac] = obj.createACBlock(initParams);

            res = vertcat(AL_res,ME_res,AC_res);
            jac = vertcat(AL_jac,ME_jac,AC_jac);
        end
        
        %% Find Jacobian Pattern
        function jac_p = getJacobianPattern(obj)
            % If Precomputed jacobian is not working well, create sparsity
            % pattern matirx and let MATLAB lsqnonlin solver to compute
            % numerical jacobians
        end
        
        %% Create Arc Length Block
        function [res,jac] = createALBlock(obj,Lf)
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
                    currL = Lf(j);
                    if currL < prevL
                        disp([i,j])
                        error('Arc Length should be strictly increasing')
                    end

                    dist = norm(obj.points(:,j) - obj.points(:,j-1));
                    if dist > 10 
                        % If points are far away from each other, straight
                        % line distance will be quite different from arc
                        % length distance.
                        cov = 0.5^2;
                    else
                        cov = 0.1^2;
                    end
                    res(cnt) = InvMahalanobis(currL - prevL - dist,cov);
                    
                    if j ~= lb+1
                        [I1,J1,V1] = sparseFormat(cnt,3+n+j-2:3+n+j-1,InvMahalanobis([-1,1],cov));
                        I(jac_cnt+1:jac_cnt+2) = I1;
                        J(jac_cnt+1:jac_cnt+2) = J1;
                        V(jac_cnt+1:jac_cnt+2) = V1;
                        jac_cnt = jac_cnt + 2;
                    else
                        [I1,J1,V1] = sparseFormat(cnt,3+n+j-1,InvMahalanobis(1,cov));
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
        function [res,jac] = createMEBlock(obj,initParams,kappa,Lf,L,nodePoints,precompX,precompY)
            n = length(kappa);
            m = size(obj.points,2)-1;
%             m = size(obj.points,2);
            blk_height = 2*m;
%             blk_heightR = 2*size(obj.points,2);

            heading = initParams(3);
            
            res = zeros(blk_height,1);
%             resR = zeros(blk_heightR,1);
            cnt = 1; 
%             cntR = 1;
            jac_cnt = 0; 
%             jac_cntR = 0;

            I = zeros(1,2*(2*sum(obj.assoc(2:end))+3*m)); J = I; V = I;
%             IR = zeros(1,2*(2*sum(obj.assoc)+3*m)); JR = IR; VR = IR;

            for i=1:n
                Idxs = find(obj.assoc == i);
%                 IdxsR = find(obj.assoc == i);
                
                xim1 = nodePoints(1,i);
                yim1 = nodePoints(2,i);                

                %% 0m Points
                if i == 1
                    jstart = 2;
                else
                    jstart = 1;
                end
                for j=jstart:length(Idxs)
                    L_ = Lf(Idxs(j));
                    cov = reshape(obj.covs(:,Idxs(j)),2,2);

                    headingN = heading + kappa(i) * L_;
                    
                    % Residual Computation
                    xij = xim1 + 1/kappa(i) * (sin(headingN) - sin(heading));
                    yij = yim1 - 1/kappa(i) * (cos(headingN) - cos(heading));

                    res(2*cnt-1:2*cnt) = InvMahalanobis([xij;yij] - obj.points(:,Idxs(j)),cov);
                    
                    % Measurement Jacobian Computation
                    init_jac = zeros(2,3);
                    init_jac(1,1) = 1; init_jac(2,2) = 1;

                    kappa_jac = zeros(2,i);
                    L_jac = zeros(2,i);
                    
                    if i == 1 % Matched with first segment
                        init_jac(1,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
                        init_jac(2,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
                        
                        % k = 1 only
                        for k=1:i
                            kappa_jac(1,k) = L_/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
                            kappa_jac(2,k) = L_/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
                            
                            L_jac(1,k) = cos(headingN);
                            L_jac(2,k) = sin(headingN);
                        end
                    else
                        init_jac(1,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
                        init_jac(2,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));
                        
                        for k=1:i
                            if k ~= i
                                kappa_jac(1,k) = precompX(i-1,3+k) + L(k)/kappa(i) * (cos(headingN) - cos(heading));
                                kappa_jac(2,k) = precompY(i-1,3+k) + L(k)/kappa(i) * (sin(headingN) - sin(heading));

                                L_jac(1,k) = precompX(i-1,3+n+k) + kappa(k)/kappa(i) * (cos(headingN) - cos(heading));
                                L_jac(2,k) = precompY(i-1,3+n+k) + kappa(k)/kappa(i) * (sin(headingN) - sin(heading));
                            else
                                kappa_jac(1,k) = L_/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
                                kappa_jac(2,k) = L_/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
                                
                                L_jac(1,k) = cos(headingN);
                                L_jac(2,k) = sin(headingN);
                            end
                        end
                    end

                    % Normalize
                    init_jac = InvMahalanobis(init_jac,cov);
                    kappa_jac = InvMahalanobis(kappa_jac,cov);
                    L_jac = InvMahalanobis(L_jac,cov);

                    % Assign
                    [I_i,J_i,V_i] = sparseFormat(2*cnt-1:2*cnt,1:3,init_jac);
                    [I_k,J_k,V_k] = sparseFormat(2*cnt-1:2*cnt,3+1:3+i,kappa_jac);
                    [I_sL,J_sL,V_sL] = sparseFormat(2*cnt-1:2*cnt,3+n+[obj.ref_idxs(1:i-1),Idxs(j)]-1,L_jac);
                    
                    I(jac_cnt+1:jac_cnt+2*(2*i+3)) = [I_i,I_k,I_sL];
                    J(jac_cnt+1:jac_cnt+2*(2*i+3)) = [J_i,J_k,J_sL];
                    V(jac_cnt+1:jac_cnt+2*(2*i+3)) = [V_i,V_k,V_sL];

                    cnt = cnt + 1;
                    jac_cnt = jac_cnt + 2*(2*i+3);
                end
                jac = sparse(I,J,V,blk_height,obj.blk_width);

                %% Rest
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
        
        %% Create Init Point Anchoring Block
        function [res,jac] = createACBlock(obj,initParams)
            X = obj.points(:,1);
%             cov = reshape(obj.covs0(:,1),2,2);
            cov = diag([1e-4,1e-4]);

            res = InvMahalanobis(X - initParams(1:2),cov);
            jac = sparse([],[],[],2,obj.blk_width);
            jac(1,1) = 1; jac(2,2) = 1;
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
            
            X0 = [initParams;kappa;Ls];
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

    end
end