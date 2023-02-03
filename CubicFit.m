classdef CubicFit < handle
% CubicFit - Incrementally perform cubic spline approximation via EKF
% 
% Very fast computational time, not continuous, too many segments
% Run this script to use as a comparison with the proposed arc spline 
% optimization framework 
%
% * CubicFit is implemented based on the model introduced at paper below
% 
% [1] G. -P. Gwon, W. -S. Hur, S. -W. Kim and S. -W. Seo, 
% "Generation of a Precise and Efficient Lane-Level Road Map for Intelligent Vehicle Systems," 
% in IEEE Transactions on Vehicular Technology, vol. 66, no. 6, pp. 4517-4533, June 2017, 
% doi: 10.1109/TVT.2016.2535210.
% 
% Here, we simplify the problem and perform estimation at 2D domain
%
% Implemented by JinHwan Jeon, 2022
%
    properties (Access = public)
        id           
        points
        covs  
        segments = {}
        Qk = diag([0.1^2,repmat([1e-2,1e-4,1e-6,1e-8],1,2)])        
    end
    methods (Access = public)
        %% Constructor
        function obj = CubicFit(points,covs,id)
            obj.points = points;
            obj.covs = covs;
            obj.id = id;
            obj.segments = {};
        end
        
        %% Optimize
        function obj = optimize(obj)
            Xkk = zeros(9,1); Xkk(2) = obj.points(1,1); Xkk(6) = obj.points(2,1);
            Pkk = diag([0,0,1e2,1e2,1e2,0,1e2,1e2,1e2]);
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
                Sk = [Sk Xkk(1)];
                point_id = [point_id, i+1];

                % Check validity of current approximation
                % If valid, continue on with current parameters
                % If not valid, end current segment and start new segment
                isValid = obj.validity(Sk,point_id,[Xkk(2:5)'; Xkk(6:9)']);
                if ~isValid
                    % Save Recent valid approximation information
                    seg = struct();
                    seg.L = prev_L;
                    seg.id = prev_id;
                    seg.paramsX = prev_paramsX;
                    seg.paramsY = prev_paramsY;
                    obj.segments = [obj.segments, {seg}];
                    
                    % Initialize new segment
                    Sk = 0;
                    point_id = i;
                    lastL = prev_L(end);
                    lastX = obj.pointProp(lastL,[prev_paramsX;prev_paramsY]);

                    Xkk = zeros(9,1);
%                     Xkk(2) = obj.points(1,i); Xkk(6) = obj.points(2,i);
                    Xkk(2) = lastX(1); Xkk(6) = lastX(2);
                    Pkk = diag([0,0,1e2,1e2,1e2,0,1e2,1e2,1e2]);

                else
                    prev_L = Sk;
                    prev_id = point_id;
                    prev_paramsX = Xkk(2:5)';
                    prev_paramsY = Xkk(6:9)';
                    i = i+1;

                    if i == size(obj.points,2)
                        seg = struct();
                        seg.L = prev_L;
                        seg.id = prev_id;
                        seg.paramsX = prev_paramsX;
                        seg.paramsY = prev_paramsY;
                        obj.segments = [obj.segments, {seg}];
                        break;
                    end
                end
            end
            
        end

        %% Visualize
        function visualize(obj)
            figure(1); hold on; grid on; axis equal;
            for i=1:length(obj.segments)
                seg = obj.segments{i};
                pc = obj.pointProp(seg.L,[seg.paramsX;seg.paramsY]);  
                
                filled_L = linspace(0,seg.L(end),1e3);
                filled_pc = obj.pointProp(filled_L,[seg.paramsX;seg.paramsY]);

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
            Xkkp1 = Xkk + [uk;zeros(8,1)];
            Fk = eye(9);
            Pkkp1 = Fk * Pkk * Fk' + obj.Qk;
        end

        %% Measurement Update
        function [Xkk,Pkk] = update(obj,Xkkp1,Pkkp1,idx)
            Zk = obj.points(:,idx); Rk = reshape(obj.covs(:,idx),2,2);
            
            Lk = Xkkp1(1); 
            C0x = Xkkp1(2); C1x = Xkkp1(3); C2x = Xkkp1(4); C3x = Xkkp1(5);
            C0y = Xkkp1(6); C1y = Xkkp1(7); C2y = Xkkp1(8); C3y = Xkkp1(9);

            Zpred = [C0x + C1x * Lk + C2x * Lk^2 + C3x * Lk^3;
                     C0y + C1y * Lk + C2y * Lk^2 + C3y * Lk^3];
            
            Hk = [C1x + 2 * C2x * Lk + 3 * C3x * Lk^2, 1, Lk, Lk^2, Lk^3, zeros(1,4);
                  C1y + 2 * C2y * Lk + 3 * C3y * Lk^2, zeros(1,4), 1, Lk, Lk^2, Lk^3];

            yk = Zk - Zpred;
            
            Sk = Hk * Pkkp1 * Hk' + Rk;
            Kk = Pkkp1 * Hk' / Sk;

            Xkk = Xkkp1 + Kk * yk;
            Pkk = (eye(9) - Kk * Hk) * Pkkp1;            
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

            if invalid_cnt >= 4 
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
            % params : [C0x C1x C2x C3x;C0y C1y C2y C3y]
            % L: [L1, ... Ln]
            PC = zeros(2,length(L));
            for i=1:length(L)
                Lmul = [1; L(i); L(i)^2; L(i)^3];
                PC(:,i) = params * Lmul;
            end
        end
       
    end
end