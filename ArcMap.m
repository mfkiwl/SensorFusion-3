classdef ArcMap < handle
    %ARCMAP - Arc Spline based lane model for SNLS optimization
    % 
    %   Detailed explanation goes here

    properties (Access = public)
        states
        lane
        segments
        arc_segments = {}
        segment_info
        lane_prob_thres
        assocL
        assocR
        dummy = struct() % Dummy variable for debugging
    end
    
    %% Public Methods
    methods (Access = public)
        %% Constructor
        function obj = ArcMap(states,lane,lane_prob_thres)
            obj.states = states;
            obj.lane = lane;
            obj.lane_prob_thres = lane_prob_thres;
            
            % Initial Segment Classification
            obj.InitSegmentClassification();
            
            % Initial Segment Parametrization
            obj.InitSegmentParametrization();

            % Data Association
            obj.associate();
        end

        %% Map Visualization (2D Segmentwise Point Map)
        function obj = visualize2DMap(obj)
            % Alter this function for debugging purpose
            % Current : Visualize initial arc segmentation results
            
%             n = length(obj.arc_segments);
%             
%             figure(1); hold on; grid on; axis equal;
        end
        
    end
    
    %% Private Methods
    methods (Access = private)
        %% Initial Segment Classification
        function obj = InitSegmentClassification(obj) 
            % Initial Segment Classification
            % Cluster "same" lane points, even for Lane Change.
            
            obj.segment_info = zeros(2,size(obj.lane.FactorValidIntvs,1));
            LeftSegNum = 1; RightSegNum = 2;
            obj.segment_info(1,1) = LeftSegNum;
            obj.segment_info(2,1) = RightSegNum;

            for i=1:length(obj.lane.LC_dirs)
                if strcmp(obj.lane.LC_dirs{i},'left')
                    tmp = LeftSegNum;
                    LeftSegNum = max([LeftSegNum,RightSegNum]) + 1;
                    RightSegNum = tmp;
                elseif strcmp(obj.lane.LC_dirs{i},'right')
                    tmp = RightSegNum;
                    RightSegNum = max([LeftSegNum, RightSegNum]) + 1;
                    LeftSegNum = tmp;
                end
                obj.segment_info(1,i+1) = LeftSegNum;
                obj.segment_info(2,i+1) = RightSegNum;            
            end

            n = size(obj.lane.FactorValidIntvs,1);
            
            obj.segments = {};

            for i=1:n+1
                obj.segments = [obj.segments {[]}];
            end

            for i=1:n
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                Left = zeros(3,ub-lb+1);
                Right = zeros(3,ub-lb+1);
                for j=lb:ub
                    R = obj.states{j}.R;
                    P = obj.states{j}.P;
                    Left(:,j-lb+1) = P + R * obj.states{j}.left(:,1);
                    Right(:,j-lb+1) = P + R * obj.states{j}.right(:,1);
                end
                
                obj.segments{obj.segment_info(1,i)} = [obj.segments{obj.segment_info(1,i)} Left];
                obj.segments{obj.segment_info(2,i)} = [obj.segments{obj.segment_info(2,i)} Right];
            end
        end
        
        %% Initial Segment Parametrization
        function obj = InitSegmentParametrization(obj)
            % For each separate segment, perform parametrization
            n = length(obj.segments);
            for i=1:n
                disp(['[Segment ',num2str(i),']'])
                initIntv = obj.segmentation(obj.segments{i},[1,length(obj.segments{i})],1);
                obj.arc_segments = [obj.arc_segments {obj.mergeSeg(obj.segments{i},initIntv)}];
            end
        end

        %% Initial segmentation
        function intvs = segmentation(obj,LP,intv,depth)
            % (1). Connect first and end intervals with a line
            % (2). Find the maximum vertical error data point
            % (3). Divide the data points of interest into two w.r.t point found at (b).
            % (4). Repeat (1)~(3) for every divided segments.
            % * If maximum error computed at (2) is lower than threshold, 
            % stop dividing current segment further and propagate backwards. 

            % Adjust this value to increase or decrease the number of line
            % interpolation. Typically, values smaller than 0.3m is
            % recommended.
            line_acc_thres = 0.2;

            %% Create Line and visualize current segmentation state 
            init_point = LP(:,intv(1));
            last_point = LP(:,intv(2));
            m = (last_point(2) - init_point(2))/(last_point(1) - init_point(1));
            n = -m*init_point(1) + init_point(2);
            % Visualize 
%             x = linspace(init_point(1),last_point(1),1e3);
%             y = m*x+n;
%             plot(x,y,'k--');
%             plot([x(1) x(end)],[y(1), y(end)],'bs')
            
            %% Find existence of intersection
            % Even if more than 1 intersection occurs, divide current segment with
            % two segments. For stability, search intersection from the middle
            X = LP(1,intv(1):intv(2)); Y = LP(2,intv(1):intv(2));
            d = (Y - m * X - n)/sqrt(m^2 + 1); % signed distance to line
            
            % If the maximum error caused by "Line Interpolation" is lower than
            % threshold, end segmentation(return current search interval). 
            % Else, divide segmentation problem using the maximum error point
            [max_val,max_idx] = max(abs(d));
%             plot(obj.LP(1,intv(1)+max_idx-1),obj.LP(2,intv(1)+max_idx-1),'cs');
%             pause(0.3);
            
            % End segmentation or divide current segment into 2 and
            % continue search (Divide and Conquer)
            if max_val < line_acc_thres
                % Current line interpolation is accurate enough
                intvs = intv;
            else
                intv1 = obj.segmentation(LP,[intv(1) intv(1)+max_idx-1],depth+1);
                intv2 = obj.segmentation(LP,[intv(1)+max_idx-1 intv(2)],depth+1);
                intvs = [intv1 intv2];
            end
        
        end
        
        %% Merge line interpolated segments for valid circular approximation
        function initSegments = mergeSeg(obj,LP,intvL)
            % Merge line intervals for valid circular approximation and
            % save data. Since covariance or weight data of lane points are
            % not available, use rmse error for determining the 'goodness'
            % of circular fitting.
            rmse_thres = 1;
            lb = 1; n = numel(intvL); cnt = 1;
            
            initSegments = {};

            while true
                
                [res,err] = CircleFit(LP(1,intvL(lb):intvL(n))',...
                                      LP(2,intvL(lb):intvL(n))',...
                                      ones(intvL(n)-intvL(lb)+1,1),...
                                      [ones(1,intvL(n)-intvL(lb)+1); ...
                                       zeros(1,intvL(n)-intvL(lb)+1); ...
                                       zeros(1,intvL(n)-intvL(lb)+1); ...
                                       ones(1,intvL(n)-intvL(lb)+1)],...
                                      0.99,false);
                if err.rmse < rmse_thres
                    disp(['Sub-Segment No. ',num2str(cnt),...
                          ' Idx ',num2str(intvL(lb)),' ~ ',num2str(intvL(n))])
                    res.bnds = [intvL(lb) intvL(n)];
                    initSegments = [initSegments res];
                    break;
                end

                max_fit_err = 0;
                % Error threshold value should be adequately large for safe
                % merging. If too large, there may not be appropriate 
                % intervals for good quality data, and if too small, short 
                % noisy intervals may be considered as a segment for this step.
                % Current sensor fusion model merge lanes almost perfectly,
                % so setting err_thres to about 1m will be fine.
                
                err_thres = 1; 
                while max_fit_err < err_thres
                    % Upper bound is found randomly between (lb+1,n) so that
                    % maximum circular fitting error is between values 3 and 5
                    ub = obj.randint(lb+1,n);                    
                    [~,err] = CircleFit(LP(1,intvL(lb):intvL(ub))',...
                                        LP(2,intvL(lb):intvL(ub))',...
                                        ones(intvL(ub)-intvL(lb)+1,1),...
                                        [ones(1,intvL(ub)-intvL(lb)+1); ...
                                         zeros(1,intvL(ub)-intvL(lb)+1); ...
                                         zeros(1,intvL(ub)-intvL(lb)+1); ...
                                         ones(1,intvL(ub)-intvL(lb)+1)],...
                                        0.99,false);
                    max_fit_err = err.emax;
                end
    
                rmse = inf;
    
                while rmse > rmse_thres
                    ub = ub - 1;
                    
                    if ub == lb
                        error('Need more line segmentation. Lower the line acc value')
                    end
                    [res,err] = CircleFit(LP(1,intvL(lb):intvL(ub))',...
                                          LP(2,intvL(lb):intvL(ub))',...
                                          ones(intvL(ub)-intvL(lb)+1,1),...
                                          [ones(1,intvL(ub)-intvL(lb)+1); ...
                                           zeros(1,intvL(ub)-intvL(lb)+1); ...
                                           zeros(1,intvL(ub)-intvL(lb)+1); ...
                                           ones(1,intvL(ub)-intvL(lb)+1)],...
                                          0.99,false);
                    rmse = err.rmse;
                    
                end
                disp(['Sub-Segment No. ',num2str(cnt),...
                      ' Idx ',num2str(intvL(lb)),' ~ ',num2str(intvL(ub))])
                res.bnds = [intvL(lb) intvL(ub)];
                initSegments = [initSegments res];
                cnt = cnt + 1;
                lb = ub;
                
                if lb == n
                    error('Need more line segmentation. Lower the line acc value')
                end
            end
        end
        
        %% Data Association
        function obj = associate(obj)
            % Data Association with current vehicle state and arc-spline
            % map
            n = size(obj.lane.FactorValidIntvs,1);
            
            obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
            obj.assocR = obj.assocL;

            for i=1:n
                % State lb and ub for each large segment
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);

                leftSegmentIdx = obj.segment_info(1,i);
                rightSegmentIdx = obj.segment_info(2,i);

                leftSegment = obj.arc_segments{leftSegmentIdx};
                rightSegment = obj.arc_segments{rightSegmentIdx};
                
                mL = length(leftSegment);
                mR = length(rightSegment);
                
                % Save arc segment node points
                leftSegPoints = []; rightSegPoints = [];
                % Left Segment
                for j=1:mL
                    seg = leftSegment{j};
                    initPoint = [seg.x + seg.R * cos(seg.th_init);
                                 seg.y + seg.R * sin(seg.th_init)];
                    leftSegPoints = [leftSegPoints initPoint];
                    if j == mL
                        lastPoint = [seg.x + seg.R * cos(seg.th_last);
                                     seg.y + seg.R * sin(seg.th_last)];
                        leftSegPoints = [leftSegPoints lastPoint];
                    end
                end

                % Right Segment
                for j=1:mR
                    seg = rightSegment{j};
                    initPoint = [seg.x + seg.R * cos(seg.th_init);
                                 seg.y + seg.R * sin(seg.th_init)];
                    rightSegPoints = [rightSegPoints initPoint];
                    if j == mL
                        lastPoint = [seg.x + seg.R * cos(seg.th_last);
                                     seg.y + seg.R * sin(seg.th_last)];
                        rightSegPoints = [rightSegPoints lastPoint];
                    end
                end

                % Perform matching for valid points
                for j=lb:ub
                    lane_idx = obj.lane.state_idxs(j);
                    % Filter Valid state idx
                    if obj.lane.prob(lane_idx,2) > obj.lane_prob_thres && obj.lane.prob(lane_idx,3) > obj.lane_prob_thres
                        R = obj.states{j}.R;
                        P = obj.states{j}.P;
                        
                        for k=1:obj.lane.prev_num
                            rel_pos = [10*(k-1);0;0];
                            pos = P + R * rel_pos;
                            
                            segIdxL = obj.match(pos,R,leftSegPoints,leftSegment);
                            segIdxR = obj.match(pos,R,rightSegPoints,rightSegment);

                            obj.assocL(j-lb+1,k) = segIdxL;
                            obj.assocR(j-lb+1,k) = segIdxR;
                        end
                        
                        
                    end
                end
            end

        end

        

    end

    %% Static Methods
    methods (Static)
        %% Get random integer with bnds
        function val = randint(val1,val2)
            % find val from [val1, val2] randomly
            val = val1-1+randi(val2-val1+1);
        end
        
        %% Point-Segment Matching
        function segIdx = match(pos,att,segPoints,segments)
            % May need to fix
            % Matching is performed in 2D manner
            n = size(segPoints,2);
            
            rpy = dcm2rpy(att);
            psi = rpy(3);
            xv = pos(1); yv = pos(2);
            
            cand = [];
            for i=1:n-1
                % Geometry based matching is performed
                % If intersection point is between the 2 node points,
                % dP will be the largest among [dP,d1,d2]
                % If dP is not the largest, intersection point is not
                % between the 2 node points, meaning that matching is not
                % appropriate.
                P1 = segPoints(:,i); P2 = segPoints(:,i+1);
                dP = (P1 - P2)' * (P1 - P2);
                
                x1 = P1(1); y1 = P1(2); 
                x2 = P2(1); y2 = P2(2);
                
                x_intersect = (x1 * (y2 - y1)/(x2 - x1) - y1 + xv * 1/tan(psi) + yv)/((y2 - y1)/(x2 - x1) + 1/tan(psi));
                y_intersect = (y2 - y1)/(x2 - x1) * (x_intersect - x1) + y1;
                P_intersect = [x_intersect;y_intersect];

                d1 = (P1 - P_intersect)' * (P1 - P_intersect);
                d2 = (P2 - P_intersect)' * (P2 - P_intersect);

                if max([dP,d1,d2]) == dP
                    % Matching Success
                    % Compute arc center to vehicle distance
                    arc_c = [segments{i}.x; segments{i}.y];
                    dist = abs(sqrt(([xv;yv] - arc_c)' * ([xv;yv] - arc_c)) - abs(1/segments{i}.kappa));
                    
                    cand = [cand [i;dist]];
                end
            end

            if ~isempty(cand)
                % If match candidate exists, pick the one with smallest
                % matching error as the correct match
                [~,idx] = min(cand(2,:));
                segIdx = cand(1,idx);
            else
                segIdx = 0;
            end
        end

    end
end