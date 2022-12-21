classdef ClothFitEKF < handle
% CLOTHFITEKF - Incrementally perform clothoid spline approximation via EKF
% 
% Just for comparison!
%
% Implemented by JinHwan Jeon, 2022
%
    properties (Access = public)
        id           
        points
        covs  
        segments = {}
        Qk = 0.3^2;       
    end
    methods (Access = public)
        %% Constructor
        function obj = ClothFitEKF(points,covs,id)
            obj.points = points;
            obj.covs = covs;
            obj.id = id;
            obj.segments = {};
        end
        
        %% Optimize
        function obj = optimize(obj)
            Xkk = zeros(6,1);
            Xkk(3) = atan2(obj.points(2,2) - obj.points(2,1),obj.points(1,2) - obj.points(1,1));
            Xkk(5) = obj.points(1,1); Xkk(6) = obj.points(2,1);
            Pkk = diag([0,1e-2,1e-6,1e-8,0,0]);
            Sk = 0;
            Xk = [Xkk(5);Xkk(6)];
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
                Xk = [Xk, [Xkk(5);Xkk(6)]];
                point_id = [point_id, i+1];

                % Check validity of current approximation
                % If valid, continue on with current parameters
                % If not valid, end current segment and start new segment
                isValid = obj.validity(Xk,point_id); % need to fix
                if ~isValid
                    % Save Recent valid approximation information
                    seg = struct();
                    seg.L = prev_L;
                    seg.id = prev_id;
                    seg.params = prev_params;    
                    seg.points = Xk;
                    obj.segments = [obj.segments, {seg}];
                    
                    % Initialize new segment
                    Sk = 0;
                    point_id = i;
                    Xkk = zeros(6,1);                    
                    Xkk(3) = atan2(obj.points(2,i+1) - obj.points(2,i),obj.points(1,i+1) - obj.points(1,i));
                    Xkk(5) = obj.points(1,i); Xkk(6) = obj.points(2,i);
                    Pkk = diag([0,1e-2,1e-6,1e-8,0,0]);
                    Xk = [Xkk(5);Xkk(6)];

                else
                    prev_L = Sk;
                    prev_id = point_id;
                    prev_params = Xkk(1:5);                    
                    i = i+1;

                    if i == size(obj.points,2)
                        seg = struct();
                        seg.L = prev_L;
                        seg.id = prev_id;
                        seg.params = prev_params;
                        seg.points = Xk;
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
                pc = seg.points;
                
%                 filled_L = linspace(0,seg.L(end),1e3);
%                 filled_pc = obj.pointProp(filled_L,seg.params);

                p_d = plot(obj.points(1,seg.id),obj.points(2,seg.id),'r.');
                p_est = plot(pc(1,:),pc(2,:),'b.');     
                
%                 p_filled = plot(filled_pc(1,:),filled_pc(2,:),'k--');
                p_est_end = plot([pc(1,1),pc(1,end)],[pc(2,1),pc(2,end)],'ms');
            end
            xlabel('Global X'); ylabel('Global Y'); title('Cubic Fitting')
            legend([p_d,p_est,p_est_end], ...
                   'Data Points','Estimated Points','Segment Nodes')
        end

    end
    methods (Access = private)  
        %% State Prediction
        function [Xkkp1,Pkkp1] = predict(obj,Xkk,Pkk,uk)
%             Xkkp1 = Xkk + [zeros(5,1);uk];
            Xkkp1 = zeros(length(Xkk),1);
            Xkkp1(1) = Xkk(1) + uk;
            Xkkp1(2:4) = Xkk(2:4);
            
            L = Xkk(1); tau0 = Xkk(2); kappa0 = Xkk(3); c = Xkk(4);
            x = Xkk(5); y = Xkk(6);
            tau = tau0 + kappa0 * L + 1/2 * c * L^2;
            kappa = kappa0 + c * L;

            Xkkp1(5) = x + cos(tau) * uk - 1/2 * kappa * sin(tau) * uk^2;
            Xkkp1(6) = y + sin(tau) * uk + 1/2 * kappa * cos(tau) * uk^2;

            [Jfx,Jfu] = obj.predJac(tau,kappa,c,L,uk);
            Pkkp1 = Jfx * Pkk * Jfx' + Jfu * obj.Qk * Jfu';
        end

        %% Measurement Update
        function [Xkk,Pkk] = update(obj,Xkkp1,Pkkp1,idx)
            Zk = obj.points(:,idx); Rk = reshape(obj.covs(:,idx),2,2);
            
            Zpred = Xkkp1(5:6);
            yk = Zk - Zpred;
            
            % Measurement Jacobian
            Hk = [0, 0, 0, 0, 1, 0;
                  0, 0, 0, 0, 0, 1];
            Sk = Hk * Pkkp1 * Hk' + Rk;
            Kk = Pkkp1 * Hk' / Sk;

            Xkk = Xkkp1 + Kk * yk;
            Pkk = (eye(6) - Kk * Hk) * Pkkp1;            
        end
        
        %% Check Validity
        function isValid = validity(obj,pc,idxs)
            % Propagate Points with current arc length and parameters
            
            invalid_cnt = 0;
            chisq_thres = chi2inv(0.99,2);
            for i=1:size(pc,2)
                Zmeas = obj.points(:,idxs(i));
                cov = reshape(obj.covs(:,idxs(i)),2,2);
                Zdiff = (Zmeas - pc(:,i));
                Ndist = Zdiff' / cov * Zdiff;
                if Ndist > chisq_thres
                    invalid_cnt = invalid_cnt + 1;
                end
            end

            if invalid_cnt >= 4
                % More than 3 points are invalid: can create new clothoid
                isValid = false;
            else
                isValid = true;
            end
        end

    end

    methods (Static)
        %% Prediction Jacobian
        function [Jfx,Jfu] = predJac(tau0,kappa0,c,L,dl)
            Jfx = eye(6);
            Jfu = zeros(6,1);

            Jfx(5,1) = - dl*sin((c*L^2)/2 + kappa0*L + tau0)*(kappa0 + L*c) - (c*dl*sin((c*L^2)/2 + kappa0*L + tau0))/2 - dl*cos((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2)*(kappa0 + L*c);
            Jfx(5,2) = - dl*sin((c*L^2)/2 + kappa0*L + tau0) - dl*cos((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2);
            Jfx(5,3) = - (dl*sin((c*L^2)/2 + kappa0*L + tau0))/2 - L*dl*sin((c*L^2)/2 + kappa0*L + tau0) - L*dl*cos((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2);
            Jfx(5,4) = - (L^2*dl*sin((c*L^2)/2 + kappa0*L + tau0))/2 - (L*dl*sin((c*L^2)/2 + kappa0*L + tau0))/2 - (L^2*dl*cos((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2))/2;

            Jfx(6,1) = dl*cos((c*L^2)/2 + kappa0*L + tau0)*(kappa0 + L*c) + (c*dl*cos((c*L^2)/2 + kappa0*L + tau0))/2 - dl*sin((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2)*(kappa0 + L*c);
            Jfx(6,2) = dl*cos((c*L^2)/2 + kappa0*L + tau0) - dl*sin((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2);
            Jfx(6,3) = (dl*cos((c*L^2)/2 + kappa0*L + tau0))/2 + L*dl*cos((c*L^2)/2 + kappa0*L + tau0) - L*dl*sin((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2);
            Jfx(6,4) = (L^2*dl*cos((c*L^2)/2 + kappa0*L + tau0))/2 + (L*dl*cos((c*L^2)/2 + kappa0*L + tau0))/2 - (L^2*dl*sin((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2))/2;
            
            Jfu(1) = 1; 
            Jfu(5) = cos((c*L^2)/2 + kappa0*L + tau0) - sin((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2);
            Jfu(6) = sin((c*L^2)/2 + kappa0*L + tau0) + cos((c*L^2)/2 + kappa0*L + tau0)*(kappa0/2 + (L*c)/2);
        end

    end
end