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
    end

    methods (Access = public)
        %% Constructor
        function obj = ArcFit2(params,points0,points,covs0,covs,id)
            obj.id = id;
            obj.params = params;
            obj.points0 = points0;
            obj.points = points;
            obj.covs0 = covs0;
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
            options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);
            X = lsqnonlin(@obj.cost_func,X0,[],[],options);

            % Save Optimization Results
            
        end

    end

    methods (Access = private)
        %% NLS Cost Function
        function [res,jac] = cost_func(obj,X0)
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
            [precompX, precompY] = obj.precomputeJac(initParams,kappa,L);

            [ME_res,ME_jac] = obj.createMEBlock(initParams,kappa,sqrtL0,sqrtL,L,nodePoints,precompX,precompY);
            [AC_res,AC_jac] = obj.createACBlock(initParams);

            res = vertcat(ME_res,AC_res);
            jac = vertact(ME_jac,AC_jac);
        end
        
        %% Create Measurement Block
        function [res,jac] = createMEBlock(obj,params)

        end

        %% Create Init Point Anchoring Block
        function [res,jac] = createACBlock(obj,initParams);
        end

        %% Data Association 0m Previewed Measurements
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
                            obj.assoc(i) = obj.assoc(idx+1);
                        else % Matched to "Prev"
                            obj.assoc(i) = obj.assoc(idx-1);
                        end
                    else
                        obj.assoc(i) = obj.assoc(idx);
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
            n = length(kappa);
            sqrtL = L.^(1/2);
            heading = initParams(3);

            precompX = zeros(n,2*n+3); precompY = precompX;            

            for i=1:n
                if i == 1
                    precompX(i,1) = 1;
                    precompX(i,3) = 1/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));
                    precompX(i,3+i) = L(i)/kappa(i) * cos(heading + kappa(i) * L(i)) - 1/kappa(i)^2 * (sin(heading + kappa(i) * L(i)) - sin(heading));
                    precompX(i,3+n+i) = cos(heading + kappa(i) * L(i)) * 2 * sqrtL(i);

                    precompY(i,2) = 1;
                    precompY(i,3) = 1/kappa(i) * (sin(heading + kappa(i) * L(i)) - sin(heading));
                    precompY(i,3+i) = L(i)/kappa(i) * sin(heading + kappa(i) * L(i)) + 1/kappa(i)^2 * (cos(heading + kappa(i) * L(i)) - cos(heading));
                    precompY(i,3+n+i) = sin(heading + kappa(i) * L(i)) * 2 * sqrtL(i);
                else                    
                    precompX(i,1) = precompX(i-1,1);
                    precompX(i,3) = precompX(i-1,3) + 1/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));
                    
                    precompY(i,2) = precompY(i-1,2);
                    precompY(i,3) = precompY(i-1,3) + 1/kappa(i) * (sin(heading + kappa(i) * L(i)) - sin(heading));

                    for j=1:i
                        if j ~= i
                            precompX(i,3+j) = precompX(i-1,3+j) + L(j)/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading));
                            precompX(i,3+n+j) = precompX(i-1,3+n+j) + kappa(j)/kappa(i) * (cos(heading + kappa(i) * L(i)) - cos(heading)) * 2 * sqrtL(j);
                        else
                            precompX(i,3+i) = L(i)/kappa(i) * cos(heading + kappa(i) * L(i)) - 1/kappa(i)^2 * (sin(heading + kappa(i) * L(i)) - sin(heading));
                            precompX(i,3+n+i) = cos(heading + kappa(i) * L(i)) * 2 * sqrtL(i);
                        end
                    end
    
                    % Yi
                end
                
                heading = heading + kappa(i) * L(i);
            end
            

        end

    end
end