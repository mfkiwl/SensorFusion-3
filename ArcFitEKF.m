classdef ArcFitEKF < handle
    properties (Access = public)
        id           
        points
        covs  
        segments = {}
        Qk = diag([1e-2,1e-2,1e-6,1e-8,0.3^2]) ;       
    end
    methods (Access = public)
        %% Constructor
        function obj = ArcFitEKF(points,covs,id)
            obj.points = points;
            obj.covs = covs;
            obj.id = id;
            obj.segments = {};
        end
        
        %% Optimize
        function obj = optimize(obj)
            Xkk = zeros(5,1); 
            Xkk(1) = obj.points(1,1); Xkk(2) = obj.points(2,1);
            Xkk(3) = atan2(obj.points(2,2) - obj.points(2,1),obj.points(1,2) - obj.points(1,1));
            Xkk(4) = 1e-2;
            Pkk = diag([0,0,1e-1,1e-2,0]);
            Sk = 0;
            point_id = 1;
            
            i = 1;
            while true
                % EKF State Prediction
                uk = norm(obj.points(:,i+1) - obj.points(:,i));
                [Xkkp1,Pkkp1] = obj.predict(Xkk,Pkk,uk);

                % EKF State Update
                [Xkk,Pkk] = obj.update(Xkkp1,Pkkp1,i+1);

                % Augment data
                Sk = [Sk Xkk(5)];
                point_id = [point_id, i+1];

                % Check validity of current approximation
                % If valid, continue on with current parameters
                % If not valid, end current segment and start new segment
                isValid = obj.validity(Sk,point_id,Xkk(1:4));
                if ~isValid
                    % Save Recent valid approximation information
                    seg = struct();
                    seg.L = prev_L;
                    seg.id = prev_id;
                    seg.params = prev_params;                    
                    obj.segments = [obj.segments, {seg}];
                    
                    % Initialize new segment
                    Sk = 0;
                    point_id = i;
                    Xkk = zeros(5,1);
                    Xkk(1) = obj.points(1,i); Xkk(2) = obj.points(2,i);
                    Xkk(3) = atan2(obj.points(2,i+1) - obj.points(2,i),obj.points(1,i+1) - obj.points(1,i));
%                     Xkk(3) = prev_params(3) + prev_params(4) * prev_L(end);
                    Xkk(4) = 1e-2;
                    Pkk = diag([0,0,1e-1,1e-2,0]);

                else
                    prev_L = Sk;
                    prev_id = point_id;
                    prev_params = Xkk(1:4);                    
                    i = i+1;

                    if i == size(obj.points,2)
                        seg = struct();
                        seg.L = prev_L;
                        seg.id = prev_id;
                        seg.params = prev_params;                        
                        obj.segments = [obj.segments, {seg}];
                        break;
                    end
                end
            end
            
        end

        %% Visualize
        function visualize(obj)
            figure(20); hold on; grid on; axis equal;
            for i=1:length(obj.segments)
                seg = obj.segments{i};
                pc = obj.pointProp(seg.L,seg.params);  
                
                filled_L = linspace(0,seg.L(end),1e3);
                filled_pc = obj.pointProp(filled_L,seg.params);

                p_d = plot(obj.points(1,seg.id),obj.points(2,seg.id),'r.');
                p_est = plot(pc(1,:),pc(2,:),'b.');     
                p_filled = plot(filled_pc(1,:),filled_pc(2,:),'k--');
                p_est_end = plot([pc(1,1),pc(1,end)],[pc(2,1),pc(2,end)],'ms');
            end
            xlabel('Global X'); ylabel('Global Y'); title('Cubic Fitting')
            legend([p_d,p_est,p_filled,p_est_end], ...
                   'Data Points','Estimated Points', ...
                   'Estimated Curve','Segment Nodes')
        end

    end
    methods (Access = private)  
        %% State Prediction
        function [Xkkp1,Pkkp1] = predict(obj,Xkk,Pkk,uk)
            Xkkp1 = Xkk + [zeros(4,1);uk];
            Fk = eye(5);
            Pkkp1 = Fk * Pkk * Fk' + obj.Qk;
        end

        %% Measurement Update
        function [Xkk,Pkk] = update(obj,Xkkp1,Pkkp1,idx)
            Zk = obj.points(:,idx); Rk = reshape(obj.covs(:,idx),2,2);
            
            x0 = Xkkp1(1); y0 = Xkkp1(2); tau0 = Xkkp1(3); kappa = Xkkp1(4);
            Lk = Xkkp1(5); 
            
            Zpred = [x0 + 1/kappa * (sin(tau0 + kappa * Lk) - sin(tau0));
                     y0 - 1/kappa * (cos(tau0 + kappa * Lk) - cos(tau0))];
            yk = Zk - Zpred;

            % Measurement Jacobian
            Hk = zeros(2,5);
            Hk(1:2,1:2) = eye(2);
            Hk(1,3) = 1/kappa * (cos(tau0 + kappa * Lk) - cos(tau0));
            Hk(1,4) = Lk/kappa * cos(tau0 + kappa * Lk) - 1/kappa^2 * (sin(tau0 + kappa * Lk) - sin(tau0));
            Hk(1,5) = cos(tau0 + kappa * Lk);

            Hk(2,3) = 1/kappa * (sin(tau0 + kappa * Lk) - sin(tau0));
            Hk(2,4) = Lk/kappa * sin(tau0 + kappa * Lk) + 1/kappa^2 * (cos(tau0 + kappa * Lk) - cos(tau0));
            Hk(2,5) = sin(tau0 + kappa * Lk);

            Sk = Hk * Pkkp1 * Hk' + Rk;
            Kk = Pkkp1 * Hk' / Sk;

            Xkk = Xkkp1 + Kk * yk;
            Pkk = (eye(5) - Kk * Hk) * Pkkp1;            
        end
        
        %% Check Validity
        function isValid = validity(obj,L,idxs,params)
            % Propagate Points with current arc length and parameters
            pc = obj.pointProp(L,params);
            invalid_cnt = 0;
            chisq_thres = chi2inv(0.99,2);
            for i=1:length(L)
                Zmeas = obj.points(:,idxs(i));
                cov = reshape(obj.covs(:,idxs(i)),2,2);
                Zdiff = (Zmeas - pc(:,i));
                Ndist = Zdiff' / cov * Zdiff;
                if Ndist > chisq_thres
                    invalid_cnt = invalid_cnt + 1;
                end
            end

            if invalid_cnt >= 3 
                % More than 4 points are invalid: can create new cubic
                % spline curve
                isValid = false;
            else
                isValid = true;
            end
        end

    end

    methods (Static)
        %% Propagate Points with given cubic parameters
        function PC = pointProp(L,params)
            % params : [x0,y0,tau0,kappa]
            % L: [L1, ... Ln]
            PC = zeros(2,length(L));
            x0 = params(1); y0 = params(2); 
            tau0 = params(3); kappa = params(4);
            for i=1:length(L)                
                PC(1,i) = x0 + 1/kappa * (sin(tau0 + kappa * L(i)) - sin(tau0));
                PC(2,i) = y0 - 1/kappa * (cos(tau0 + kappa * L(i)) - cos(tau0));
            end
        end
       
    end
end