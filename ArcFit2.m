classdef ArcFit2 < handle
    properties(Access = public)
        id
        params        
        points0
        points
        covs0
        covs                
        assoc0
        assoc
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
            obj.associate0();
            obj.associate();
        end
        
        %% Retrieve params 
        function params = getParams(obj)
            params = obj.params;
        end
        
        %% Optimizer
        function obj = optimize(obj)
            X0 = obj.createInitVal();
            obj.blk_width = length(X0);
            jac_p = sparse([0,0,0,0,1,0]);
            options = optimoptions('lsqnonlin', ...                                                                      
                                   'MaxFunctionEvaluations',2000, ...
                                   'MaxIterations',2000, ...
                                   'Display','iter-detailed', ...
                                   'SpecifyObjectiveGradient',true);
            disp(['[Performing Optimization for Segment ID: ',num2str(obj.id),']']) 
            X = lsqnonlin(@obj.cost_func,X0,[],[],options);

            % Save Optimization Results
            n = length(obj.params.kappa);
            obj.params.x0 = X(1);
            obj.params.y0 = X(2);
            obj.params.tau0 = X(3);
            obj.params.kappa = X(3+1:3+n)';
            
            sqrtL = X(3+n+obj.ref_idxs)';
            obj.params.L = sqrtL.^2;

        end

    end

    methods (Access = private)
        %% NLS Cost Function
        function res = cost_func(obj,X0)
            n = length(obj.params.kappa);
            initParams = X0(1:3);
            kappa = X0(3+1:3+n)';
            sqrtL0 = X0(3+n+1:3+n+size(obj.points0,2))';
            sqrtL = X0(3+n+size(obj.points0,2)+1:end)';
            
            L = zeros(1,n);
            for i=1:size(obj.params.bnds,1)
                L(i) = (sqrtL0(obj.params.bnds(i,2)))^2;
            end

            nodePoints = obj.propNode(initParams,kappa,L);
            % Pre compute Jacobians
%             [precompX, precompY] = obj.precomputeJac(initParams,kappa,L);

            ME_res = obj.createMEBlock(initParams,kappa,sqrtL0,sqrtL,L,nodePoints);
            AC_res = obj.createACBlock(initParams);

            res = vertcat(ME_res,AC_res);
%             jac = vertcat(ME_jac,AC_jac);
        end
        
        %% Find Jacobian Pattern
        function jac_p = getJacobianPattern(obj)
            
        end

        %% Create Measurement Block
        function res = createMEBlock(obj,initParams,kappa,sqrtL0,sqrtL,L,nodePoints)
            n = length(kappa);
%             m0 = size(obj.points0,2);
%             m = size(obj.points,2);
            blk_height0 = 2*size(obj.points0,2);
            blk_heightR = 2*size(obj.points,2);

            heading = initParams(3);
            
            res0 = zeros(blk_height0,1);
            resR = zeros(blk_heightR,1);
            cnt0 = 1; cntR = 1;
%             jac_cnt0 = 0; jac_cntR = 0;

%             I0 = zeros(1,2*(2*sum(obj.assoc0)+3*m0)); J0 = I0; V0 = I0;
%             IR = zeros(1,2*(2*sum(obj.assoc)+3*m)); JR = IR; VR = IR;

            for i=1:n
                Idxs0 = find(obj.assoc0 == i);
                IdxsR = find(obj.assoc == i);
                
                xim1 = nodePoints(1,i);
                yim1 = nodePoints(2,i);                

                %% 0m Points
                for j=1:length(Idxs0)
                    sqrtL0_ = sqrtL0(Idxs0(j));
                    cov = reshape(obj.covs0(:,Idxs0(j)),2,2);

                    headingN = heading + kappa(i) * sqrtL0_^2;
                    
                    % Residual Computation
                    xij = xim1 + 1/kappa(i) * (sin(headingN) - sin(heading));
                    yij = yim1 - 1/kappa(i) * (cos(headingN) - cos(heading));

                    res0(2*cnt0-1:2*cnt0) = InvMahalanobis([xij;yij] - obj.points0(:,Idxs0(j)),cov);
                    
%                     % Measurement Jacobian Computation
%                     init_jac0 = zeros(2,3);
%                     init_jac0(1,1) = 1; init_jac0(2,2) = 1;
% 
%                     kappa_jac0 = zeros(2,i);
%                     sqrtL_jac0 = zeros(2,i);
%                     
%                     if i == 1 % Matched with first segment
%                         init_jac0(1,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
%                         init_jac0(2,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
%                         
%                         % k = 1 only
%                         for k=1:i
%                             kappa_jac0(1,k) = sqrtL0_^2/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                             kappa_jac0(2,k) = sqrtL0_^2/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                             
%                             sqrtL_jac0(1,k) = cos(headingN) * 2 * sqrtL0_;
%                             sqrtL_jac0(2,k) = sin(headingN) * 2 * sqrtL0_;
%                         end
%                     else
%                         init_jac0(1,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
%                         init_jac0(2,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));
%                         
%                         for k=1:i
%                             if k ~= i
%                                 kappa_jac0(1,k) = precompX(i-1,3+k) + L(k)/kappa(i) * (cos(headingN) - cos(heading));
%                                 kappa_jac0(2,k) = precompY(i-1,3+k) + L(k)/kappa(i) * (sin(headingN) - sin(heading));
% 
%                                 sqrtL_jac0(1,k) = precompX(i-1,3+n+k) + kappa(k)/kappa(i) * (cos(headingN) - cos(heading)) * 2 * sqrtL0_;
%                                 sqrtL_jac0(2,k) = precompY(i-1,3+n+k) + kappa(k)/kappa(i) * (sin(headingN) - sin(heading)) * 2 * sqrtL0_;
%                             else
%                                 kappa_jac0(1,k) = sqrtL0_^2/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
%                                 kappa_jac0(2,k) = sqrtL0_^2/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
%                                 
%                                 sqrtL_jac0(1,k) = cos(headingN) * 2 * sqrtL0_;
%                                 sqrtL_jac0(2,k) = sin(headingN) * 2 * sqrtL0_;
%                             end
%                         end
%                     end
% 
%                     % Normalize
%                     init_jac0 = InvMahalanobis(init_jac0,cov);
%                     kappa_jac0 = InvMahalanobis(kappa_jac0,cov);
%                     sqrtL_jac0 = InvMahalanobis(sqrtL_jac0,cov);
% 
%                     % Assign
%                     [I0_i,J0_i,V0_i] = sparseFormat(2*cnt0-1:2*cnt0,1:3,init_jac0);
%                     [I0_k,J0_k,V0_k] = sparseFormat(2*cnt0-1:2*cnt0,3+1:3+i,kappa_jac0);
%                     [I0_sL,J0_sL,V0_sL] = sparseFormat(2*cnt0-1:2*cnt0,3+n+[obj.ref_idxs(1:i-1),Idxs0(j)],sqrtL_jac0);
%                     
%                     I0(jac_cnt0+1:jac_cnt0+2*(2*i+3)) = [I0_i,I0_k,I0_sL];
%                     J0(jac_cnt0+1:jac_cnt0+2*(2*i+3)) = [J0_i,J0_k,J0_sL];
%                     V0(jac_cnt0+1:jac_cnt0+2*(2*i+3)) = [V0_i,V0_k,V0_sL];

                    cnt0 = cnt0 + 1;
%                     jac_cnt0 = jac_cnt0 + 2*(2*i+3);
                end

                %% Rest
                for j=1:length(IdxsR)
                    sqrtL_ = sqrtL(IdxsR(j));
                    cov = reshape(obj.covs(:,IdxsR(j)),2,2);

                    headingN = heading + kappa(i) * sqrtL_^2;
                    
                    % Residual Computation
                    xij = xim1 + 1/kappa(i) * (sin(headingN) - sin(heading));
                    yij = yim1 - 1/kappa(i) * (cos(headingN) - cos(heading));

                    resR(2*cntR-1:2*cntR) = InvMahalanobis([xij;yij] - obj.points(:,IdxsR(j)),cov);
                    
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

                    cntR = cntR + 1;
%                     jac_cntR = jac_cntR + 2*(2*i+3);
                end

                heading = heading + kappa(i) * L(i);
            end

%             jac0 = sparse(I0,J0,V0,blk_height0,obj.blk_width);
%             jacR = sparse(IR,JR,VR,blk_heightR,obj.blk_width);
            
            res = vertcat(res0,resR);
%             jac = vertcat(jac0,jacR);

        end
        
        %% Create Init Point Anchoring Block
        function res = createACBlock(obj,initParams)
            X = obj.points0(:,1);
%             cov = reshape(obj.covs0(:,1),2,2);
            cov = diag([1e-4,1e-4]);

            res = InvMahalanobis(X - initParams(1:2),cov);
%             jac = sparse([],[],[],2,obj.blk_width);
%             jac(1,1) = 1; jac(2,2) = 1;
%             jac = InvMahalanobis(jac,cov);
        end

        %% Data Association 0m Previewed Measurements
        function obj = associate0(obj)
            obj.assoc0 = zeros(1,size(obj.points0,2));
            % If state idx is used, then association is fixed throughout
            % the whole optimization process
            bnds = obj.params.bnds;
            for i=1:size(bnds,1)
                lb = bnds(i,1); ub = bnds(i,2);
                if i~=1
                    obj.assoc0(lb+1:ub) = i;
                else
                    obj.assoc0(lb:ub) = i;
                end
            end

            % Create reference idxs for (end) node points
            obj.ref_idxs = obj.params.bnds(:,2)';
        end

        %% Data Association Remaining Measurements
        function obj = associate(obj)
            %ASSOCIATE Data Association 
            % Matches which point belongs to which segment
            obj.assoc = zeros(1,size(obj.points,2));            
            for i=1:size(obj.points,2)
                Px = obj.points(1,i); Py = obj.points(2,i);
                P_refx = obj.points0(1,:); P_refy = obj.points0(2,:);
                D = sqrt((P_refx - Px).^2 + (P_refy - Py).^2);
                [~,idx] = min(D);
                
                if idx ~= 1 && idx ~= size(obj.points0,2)
                    % Check if idx belongs to any of the boundary index
                    [r,~] = find(obj.params.bnds == idx);
                    
                    % If nearest 0m point is matched, then determine
                    % association using nearby 0m points
                    if ~isempty(r)
                        prev = obj.points0(idx-1);
                        next = obj.points0(idx+1);
                        
                        d1 = norm([Px;Py] - prev);
                        d2 = norm([Px;Py] - next);

                        if d1 >= d2 % Matched to "Next"
                            obj.assoc(i) = obj.assoc0(idx+1);
                        else % Matched to "Prev"
                            obj.assoc(i) = obj.assoc0(idx-1);
                        end
                    else
                        obj.assoc(i) = obj.assoc0(idx);
                    end
                else
                    obj.assoc(i) = obj.assoc0(idx);
                end
            end                
        end
        
        %% Create Initial Value
        function X0 = createInitVal(obj)
            initParams = [obj.params.x0; obj.params.y0; obj.params.tau0];
            kappa = obj.params.kappa';
            nodeSqrtL = (obj.params.L).^(1/2);

            % L: 0m L + Remaining L
            % Use Association information to assign initial value
            sqrtL0 = zeros(size(obj.points0,2),1);
            sqrtL = zeros(size(obj.points,2),1);
            
            for i=1:size(obj.points0,2)
                SegIdx = obj.assoc0(i);
                if i == obj.params.bnds(SegIdx,2)
                    % Initial value for node points
                    % sqrt(L)
                    sqrtL0(i) = nodeSqrtL(SegIdx);
                else
                    % Initial value for other points
                    % sqrt(L/2)
                    sqrtL0(i) = sqrt(nodeSqrtL(SegIdx)^2/2);
                end
            end
            
            for i=1:size(obj.points,2)
                SegIdx = obj.assoc(i);
                sqrtL(i) = sqrt(nodeSqrtL(SegIdx)^2/2);
            end
            
            X0 = [initParams;kappa;sqrtL0;sqrtL];
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
            sqrtL = L.^(1/2);
            heading = initParams(3);

            precompX = zeros(n,2*n+3); precompY = precompX;            

            for i=1:n
                headingN = heading + kappa(i) * L(i);
                if i == 1
                    precompX(i,1) = 1;
                    precompX(i,3) = 1/kappa(i) * (cos(headingN) - cos(heading));
                    precompX(i,3+i) = L(i)/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
                    precompX(i,3+n+i) = cos(headingN) * 2 * sqrtL(i);

                    precompY(i,2) = 1;
                    precompY(i,3) = 1/kappa(i) * (sin(headingN) - sin(heading));
                    precompY(i,3+i) = L(i)/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
                    precompY(i,3+n+i) = sin(headingN) * 2 * sqrtL(i);
                else                    
                    precompX(i,1) = precompX(i-1,1);
                    precompX(i,3) = precompX(i-1,3) + 1/kappa(i) * (cos(headingN) - cos(heading));
                    
                    precompY(i,2) = precompY(i-1,2);
                    precompY(i,3) = precompY(i-1,3) + 1/kappa(i) * (sin(headingN) - sin(heading));

                    for j=1:i
                        if j ~= i
                            precompX(i,3+j) = precompX(i-1,3+j) + L(j)/kappa(i) * (cos(headingN) - cos(heading));
                            precompX(i,3+n+j) = precompX(i-1,3+n+j) + kappa(j)/kappa(i) * (cos(headingN) - cos(heading)) * 2 * sqrtL(j);
                            
                            precompY(i,3+j) = precompY(i-1,3+j) + L(j)/kappa(i) * (sin(headingN) - sin(heading));
                            precompY(i,3+n+j) = precompY(i-1,3+n+j) + kappa(j)/kappa(i) * (sin(headingN) - sin(heading)) * 2 * sqrtL(j);
                        else
                            precompX(i,3+i) = L(i)/kappa(i) * cos(headingN) - 1/kappa(i)^2 * (sin(headingN) - sin(heading));
                            precompX(i,3+n+i) = cos(headingN) * 2 * sqrtL(i);

                            precompY(i,3+i) = L(i)/kappa(i) * sin(headingN) + 1/kappa(i)^2 * (cos(headingN) - cos(heading));
                            precompY(i,3+n+i) = sin(headingN) * 2 * sqrtL(i);                        
                        end
                    end                        
                end
                
                heading = headingN;
            end
        end

    end
end