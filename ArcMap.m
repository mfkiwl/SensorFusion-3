classdef ArcMap < handle
    %ARCMAP - Arc Spline based lane model for SNLS optimization
    % 
    %   Detailed explanation goes here

    properties (Access = public)
        states
        lane
        segments
        arc_segments = {}
        arc_nodes = {}
        segment_info
        lane_prob_thres
        assocL
        assocR
        subseg_cnt = []
        validity
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
            obj.associate();% --> initial data association is done using init segmentation information
        end
        
        %% Map Update
        function obj = update(obj,states,arc_delta_params)
            % Update variables after iterative optimization step
            % Input variable format
            % states: cell-struct --> copy of state variable
            % arc_params: cell-struct --> processed 
            obj.states = states;
            
            if length(obj.arc_segments) ~= length(arc_delta_params)
                error('Length of processed arc parameters are inappropriate')
            end
            
            % Re-initialize variables
            % Need to be careful if number of arc segments is changed
            for i=1:length(obj.arc_segments)                
                obj.arc_segments{i}.kappa = obj.arc_segments{i}.kappa + arc_delta_params{i}.kappa;
                obj.arc_segments{i}.L = obj.arc_segments{i}.L + arc_delta_params{i}.L;
                obj.arc_segments{i}.x0 = obj.arc_segments{i}.x0 + arc_delta_params{i}.x0;
                obj.arc_segments{i}.y0 = obj.arc_segments{i}.y0 + arc_delta_params{i}.y0;
                obj.arc_segments{i}.tau0 = obj.arc_segments{i}.tau0 + arc_delta_params{i}.tau0;
                obj.subseg_cnt(i) = length(arc_delta_params{i}.kappa);
            end

            % Perform data association with updated variables
%             obj.associate();
        end

        %% Map Visualization (2D Segmentwise Point Map)
        function obj = visualize2DMap(obj)
            % Alter this function for debugging purpose
            % Current : Visualize initial arc segmentation results

            figure(1); hold on; grid on; axis equal;
            
            n = length(obj.states);
            
            for i=1:n
                P = obj.states{i}.P;
                plot(P(1),P(2),'r.');
            end

            n = length(obj.arc_segments);
            for i=1:n
                m = length(obj.arc_segments{i}.kappa);
                seg = obj.arc_segments{i};
                heading = seg.tau0;
                segPoints = [seg.x0;
                             seg.y0];
                for j=1:m
                    kappa = seg.kappa(j); L = seg.L(j);
                    headingPrev = heading;
                    heading = heading + kappa * L;
                    headingCurr = heading;

                    heading_ = linspace(headingPrev,headingCurr,1e3);
                    segPoints = [segPoints segPoints(:,end) + 1/kappa * [sin(heading_) - sin(headingPrev);
                                                                        -cos(heading_) + cos(headingPrev)]];
                end
                plot(segPoints(1,:),segPoints(2,:),'k-');
            end
           
        end
        
        %% Validate Current Map after convergence
        function obj = validate(obj,phase)
            % For each segment, find the most invalid subsegment
            n = size(obj.lane.FactorValidIntvs,1);
            
            % Validity: stores number of invalid measurements for each
            % sub-segment for every segment
            obj.validity = {};
            for i=1:length(obj.arc_segments)
                obj.validity = [obj.validity {zeros(1,length(obj.arc_segments{i}.kappa))}];
            end

            for i=1:n
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                
                leftSegIdx = obj.segment_info(1,i);
                rightSegIdx = obj.segment_info(2,i);

                for j=lb:ub
                    if phase == 2
                        prev_max = 1;
                    elseif phase == 3
                        prev_max = obj.lane.prev_num;
                    end
                    
                    for k=1:prev_max
                        leftsubSegIdx = obj.assocL(j,k);
                        rightsubSegIdx = obj.assocR(j,k);

                        if leftsubSegIdx > 0
                            valid = obj.isValid(leftSegIdx,leftsubSegIdx,j,k,'left');
                            if ~valid
                                obj.validity{leftSegIdx}(leftsubSegIdx) = obj.validity{leftSegIdx}(leftsubSegIdx) + 1;
                            end
                        end

                        if rightsubSegIdx > 0
                            valid = obj.isValid(rightSegIdx,rightsubSegIdx,j,k,'right');
                            if ~valid
                                obj.validity{rightSegIdx}(rightsubSegIdx) = obj.validity{rightSegIdx}(rightsubSegIdx) + 1;
                            end
                        end
                    end
                end
            end

            % Add new segments for the most "invalid" sub-segment
            for i=1:length(obj.validity)
                [~,badsubSegIdx] = max(obj.validity{i});
                obj.replicate(i,badsubSegIdx);
            end

            obj.associate();
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
                initSegments = obj.mergeSeg(obj.segments{i},initIntv);
                initParams = struct();
                initParams.kappa = [];
                initParams.L = [];
                initParams.bnds = [];
                
                initStateIdx = obj.lane.FactorValidIntvs(obj.findColIdx(obj.segment_info,i),1);
                initParams.x0 = obj.segments{i}(1,1);
                initParams.y0 = obj.segments{i}(2,1);
                rpy = dcm2rpy(obj.states{initStateIdx}.R);
                initParams.tau0 = rpy(3); 
                heading = initParams.tau0;

                for j=1:length(initSegments)
                    initParams.bnds = [initParams.bnds; initSegments{j}.bnds];
                    initParams.kappa = [initParams.kappa initSegments{j}.kappa];
                    initParams.L = [initParams.L initSegments{j}.L];
                    if j == 1
                        initParams.xc = initParams.x0 - 1/initParams.kappa(end) * sin(heading);
                        initParams.yc = initParams.y0 + 1/initParams.kappa(end) * cos(heading);
                    else
                        initParams.xc = [initParams.xc initParams.xc(end) + (1/initParams.kappa(end-1) - 1/initParams.kappa(end)) * sin(heading)];
                        initParams.yc = [initParams.yc initParams.yc(end) - (1/initParams.kappa(end-1) - 1/initParams.kappa(end)) * cos(heading)];
                    end
                    heading = heading + initParams.kappa(end) * initParams.L(end);
                end
                obj.arc_segments = [obj.arc_segments {initParams}];
                obj.subseg_cnt = [obj.subseg_cnt length(initParams.kappa)];
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
            % parameters. After each step in optimization, vehicle states
            % and arc parameters should be updated before re-association
            %
            % Create arc segment node points before data association
            n = length(obj.arc_segments);
            for i=1:n
                m = length(obj.arc_segments{i}.kappa);
                seg = obj.arc_segments{i};
                heading = seg.tau0;
                segPoints = [seg.x0;
                             seg.y0];
                for j=1:m                    
                    kappa = seg.kappa(j); L = seg.L(j);
                    headingPrev = heading;
                    heading = heading + kappa * L;
                    headingCurr = heading;
                    segPoints = [segPoints segPoints(:,end) + 1/kappa * [sin(headingCurr) - sin(headingPrev);
                                                                        -cos(headingCurr) + cos(headingPrev)]];
                end
                obj.arc_nodes = [obj.arc_nodes {segPoints}];
            end

            % Data Association
            n = size(obj.lane.FactorValidIntvs,1);
            
            obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
            obj.assocR = obj.assocL;

            for i=1:n
                % State lb and ub for each large segment
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);

                leftSegmentIdx = obj.segment_info(1,i);
                rightSegmentIdx = obj.segment_info(2,i);

                leftSegPoints = obj.arc_nodes{leftSegmentIdx};
                rightSegPoints = obj.arc_nodes{rightSegmentIdx};
                leftSegment = obj.arc_segments{leftSegmentIdx};
                rightSegment = obj.arc_segments{rightSegmentIdx};

                % Perform matching for valid points
                for j=lb:ub
                    lane_idx = obj.lane.state_idxs(j);
                    
                    R = obj.states{j}.R;
                    P = obj.states{j}.P;
                    % Filter Valid state idx
                    if obj.lane.prob(lane_idx,2) > obj.lane_prob_thres && obj.lane.prob(lane_idx,3) > obj.lane_prob_thres
                        for k=1:obj.lane.prev_num
                            rel_pos = [10*(k-1);0;0];
                            pos = P + R * rel_pos;
                            
                            segIdxL = obj.match(pos,R,leftSegPoints,leftSegment);
                            segIdxR = obj.match(pos,R,rightSegPoints,rightSegment);

                            obj.assocL(j,k) = segIdxL;
                            obj.assocR(j,k) = segIdxR;
                        end
                    else % For low reliability lane measurement, use only 0m preview measurement                        
                        segIdxL = obj.match(P,R,leftSegPoints,leftSegment);
                        segIdxR = obj.match(P,R,rightSegPoints,rightSegment);

                        obj.assocL(j,1) = segIdxL;
                        obj.assocR(j,1) = segIdxR;
                    end
                end
            end
        end
        
        %% Check validity for each specific cases
        function flag = isValid(obj,segIdx,subsegIdx,state_idx,prev_idx,dir)
            seg = obj.arc_segments{segIdx};
            kappa = seg.kappa(1:subsegIdx); L = seg.L(1:subsegIdx);
            x0 = seg.x0; y0 = seg.y0; tau0 = seg.tau0;

            [xc,yc] = obj.getCenter(x0,y0,tau0,kappa,L);
            R = obj.states{state_idx}.R; 
            P = obj.states{state_idx}.P;
            
            P_prop = P + R * [10*(prev_idx-1);0;0];
            rpy = dcm2rpy(R); psi = rpy(3);
            R2d = [cos(psi) -sin(psi); sin(psi) cos(psi)];
            
            Xc_b = R2d' * ([xc;yc] - P_prop(1:2));
            xc_b = Xc_b(1); yc_b = Xc_b(2);


            lane_idx = obj.lane.state_idxs(state_idx);
            if strcmp(dir,'left')
                dij = obj.lane.ly(lane_idx,prev_idx);
                cov = obj.lane.lystd(lane_idx,prev_idx)^2;
            elseif strcmp(dir,'right')
                dij = obj.lane.ry(lane_idx,prev_idx);
                cov = obj.lane.rystd(lane_idx,prev_idx)^2;
            end

            if kappa(subsegIdx) > 0
                Z_pred = yc_b - sqrt((1/kappa(subsegIdx))^2 - xc_b^2);
            else
                Z_pred = yc_b + sqrt((1/kappa(subsegIdx))^2 - xc_b^2);
            end

            thres = chi2inv(0.99,2);

            chisq = (Z_pred-dij)' / cov * (Z_pred-dij); 
            if chisq > thres
                flag = false;
            else
                flag = true;
            end            
        end
        
        %% Replicate invalid segment for new optimization
        function obj = replicate(obj,segIdx,subsegIdx)
            seg = obj.arc_segments{segIdx};
            kappa = seg.kappa(subsegIdx); L = seg.L(subsegIdx);

            % Create replica by halving the original arc length
            if subsegIdx == 1
                seg.kappa = [kappa seg.kappa];
                seg.L(subsegIdx) = 1/2 * seg.L(subsegIdx);
                seg.L = [L/2 seg.L];
            elseif subsegIdx == length(seg.kappa)
                seg.kappa = [seg.kappa kappa];
                seg.L(subsegIdx) = 1/2 * seg.L(subsegIdx);
                seg.L = [seg.L L/2];
            else
                kappaF = seg.kappa(1:subsegIdx-1);
                kappaB = seg.kappa(subsegIdx+1:end);
                LF = seg.L(1:subsegIdx-1);
                LB = seg.L(subsegIdx+1:end);
                seg.kappa = [kappaF kappa kappa kappaB];
                seg.L = [LF 1/2*L 1/2 * L LB];
            end

            % Update segment info
            obj.arc_segments{segIdx} = seg;
            obj.subseg_cnt(segIdx) = obj.subseg_cnt(segIdx) + 1;
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
            
            rpy = dcm2rpy(att); psi = rpy(3);
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
                [~,idx] = max([dP,d1,d2]);
                if idx == 1
                    % Matching Success
                    % Compute arc center to vehicle distance

                    arc_c = [segments.xc(i); segments.yc(i)];
                    dist = abs(sqrt(([xv;yv] - arc_c)' * ([xv;yv] - arc_c)) - abs(1/segments.kappa(i)));
                    
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
        
        %% Compute center points of arcs, given parameters
        function [xc,yc] = getCenter(x0,y0,tau0,kappa,L)
            heading = tau0;
            xc = []; yc = [];
            for j=1:length(kappa)
                if j == 1
                    xc = x0 - 1/kappa(j) * sin(heading);
                    yc = y0 + 1/kappa(j) * cos(heading);
                else
                    xc = [xc xc(end) + (1/kappa(j-1) - 1/kappa(j)) * sin(heading)];
                    yc = [yc yc(end) - (1/kappa(j-1) - 1/kappa(j)) * cos(heading)];
                end
                heading = heading + kappa(j) * L(j);
            end
        end
        
        %% Find minimum column index for first appearance
        function idx = findColIdx(arr,num)
            % Find the first column in arr where num first appears
            [~,c] = find(arr==num);
            idx = min(c);
        end

    end
end