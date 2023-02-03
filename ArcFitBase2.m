classdef ArcFitBase2 < handle
%ARCFITBASE ArcFitBase2 performs rough optimization to get 'accurate enough'
% initial value for ArcFit. 
% 'Accurate Enough' means that at least the optimized curve passes through
% node points with acceptable accuracy precision.
% Incremental Optimization
    properties
        points % Node Points + Each Segments Center point
        params % Parameters
        m % Number of extra points per segment for stable arc fit
    end

    methods (Access = public)
        %% Constructor
        function obj = ArcFitBase2(points,params)
            obj.points = points;
            obj.params = params;
        end
        
        %% Optimize
        function obj = optimize(obj)
            disp('Running ArcFitBase for stabilizing')
            for i=1:length(obj.params.kappa)
                disp(['Current number of sub-segments: ',num2str(i)])
                X0 = [obj.params.x0; obj.params.y0; obj.params.tau0; 1./obj.params.kappa(1:i)'; obj.params.L(1:i)'];
                jac_pattern = obj.getJacPattern(i);
                options = optimoptions('lsqnonlin', ...
                                       'UseParallel',true, ...
                                       'Display','iter-detailed', ...
                                       'MaxFunctionEvaluations',inf, ...    
                                       'MaxIterations',1e4, ...
                                       'FiniteDifferenceType','forward', ...
                                       'JacobPattern',jac_pattern);
                X = lsqnonlin(@obj.cost_func,X0,[],[],options);
                obj.params.x0 = X(1);
                obj.params.y0 = X(2);
                obj.params.tau0 = X(3);
                obj.params.kappa(1:i) = 1./X(3+1:3+i)';
                obj.params.L(1:i) = X(3+i+1:end)';
            end

%             X0 = [obj.params.x0; obj.params.y0; obj.params.tau0; 1./obj.params.kappa'; obj.params.L'];
%             jac_pattern = obj.getJacPattern();
%             options = optimoptions('lsqnonlin', ...
%                                    'UseParallel',true, ...
%                                    'Display','iter-detailed', ...
%                                    'MaxFunctionEvaluations',inf, ...    
%                                    'MaxIterations',1e4, ...
%                                    'FiniteDifferenceType','forward', ...
%                                    'JacobPattern',jac_pattern);
%             X = lsqnonlin(@obj.cost_func,X0,[],[],options);
%             n = length(obj.params.kappa);
%             obj.params.x0 = X(1);
%             obj.params.y0 = X(2);
%             obj.params.tau0 = X(3);
%             obj.params.kappa = 1./X(3+1:3+n)';
%             obj.params.L = X(3+n+1:end)';
        end
        
        %% Visualize
        function obj = visualize(obj)
            n = length(obj.params.kappa);
            nodePoints = obj.propNode([obj.params.x0,obj.params.y0,obj.params.tau0], ...
                                      obj.params.kappa,obj.params.L);
            figure(1); grid on; hold on; axis equal;
            plot(obj.points(1,1:n+1),obj.points(2,1:n+1),'ms');
            plot(nodePoints(1,:),nodePoints(2,:),'cs');
            plot(obj.points(1,n+1+1:end),obj.points(2,n+1+1:end),'gx');
            
            heading = obj.params.tau0;
            SegPoints = [];
            for i=1:n
                kappa = obj.params.kappa(i); L = obj.params.L(i);
                headingPrev = heading;
                heading = heading + kappa * L;
                headingCurr = heading;

                heading_ = linspace(headingPrev,headingCurr,1e3);
                addedSegPoints = nodePoints(:,i) + 1/kappa * [sin(heading_) - sin(headingPrev);
                                                              -cos(heading_) + cos(headingPrev)];
                SegPoints = [SegPoints, addedSegPoints];
            end
            plot(SegPoints(1,:),SegPoints(2,:),'k--');

        end

    end

    methods (Access = private)
        %% NLS Cost Function
        function res = cost_func(obj,x0)
            L = length(x0); n = (L-3)/2;
%             n = length(obj.params.kappa);
            initParams = x0(1:3)';
            kappa = 1./x0(3+1:3+n)';
            L = x0(3+n+1:end)';
            Init_res = obj.CreateInitBlock(initParams);
            AM_res = obj.CreateAMBlock(initParams,kappa,L);
            ME_res = obj.CreateMEBlock(initParams,kappa,L);
            res = vertcat(Init_res,AM_res,ME_res);
        end

        %% Jacobian Pattern for faster jacbian matrix computation
        function jac_pattern = getJacPattern(obj,numSubSeg)
%             n = length(obj.params.kappa);
            blk_width = 3 + 2 * numSubSeg;
            blk_height = 2 + 2*numSubSeg + size(obj.points,2) - 1;
            jac_pattern = zeros(blk_height,blk_width);
            bnds = obj.params.bnds;
            % Init Block
            jac_pattern(1:2,1:2) = eye(2);

            % AM Block
            for i=1:numSubSeg
                jac_pattern(2+2*i-1:2+2*i,1:3+i) = ones(2,3+i);
            end

            % ME Block
            cnt = 0;
            for i=1:numSubSeg
                lb = bnds(i,1)+1; ub = bnds(i,2);
                idxs = [1:3+i,3+numSubSeg+1:3+numSubSeg+i];
                jac_pattern(2+2*numSubSeg+cnt+1:2+2*numSubSeg+cnt+ub-lb+1,idxs) = ones(ub-lb+1,3+2*i);
                cnt = cnt + ub - lb + 1;
            end
        end

        %% Create Init Block
        function res = CreateInitBlock(obj,initParams)
            cov = 1e-4 * eye(2);
            refPoint = obj.points(:,1);
            res = InvMahalanobis([initParams(1);initParams(2)] - refPoint,cov);
        end

        %% Create Anchoring Measurement Block
        function res = CreateAMBlock(obj,initParams,kappa,L)
            nodePoints = obj.propNode(initParams,kappa,L);
            n = length(kappa);
            res = zeros(2*n,1);
            cov = 1e-4 * eye(2);
            bnds = obj.params.bnds;
            for i=1:n                
                refPoint = obj.points(:,bnds(i,2));
                res(2*i-1:2*i) = InvMahalanobis(nodePoints(:,i+1) - refPoint,cov);
            end
        end

        %% Create Measurement Block
        function res = CreateMEBlock(obj,initParams,kappa,L)
            centerPoints = obj.propCenter(initParams,kappa,L);
            n = length(kappa);
            res = zeros(size(obj.points,2)-1,1);
            cov = 1e-3;
            bnds = obj.params.bnds;
            for i=1:n
                Xc = centerPoints(:,i);
                lb = bnds(i,1)+1; ub = bnds(i,2);
                for j=lb:ub
                    refPoint = obj.points(:,j);
                    res(j-1) = InvMahalanobis(norm(Xc - refPoint) - 1/abs(kappa(i)),cov);
                end
            end
        end
        
        %% Create Arc Length Anchor Block
        function res = CreateALBlock(obj,L)
            n = length(obj.params.L);
            res = zeros(n,1);
            for i=1:n
                res(i) = L(i) - obj.params.L(i);
            end
        end

    end

    methods (Static)
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

    end
end