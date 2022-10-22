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
        arc_centers = {}
        segment_info
        lane_prob_thres
        assocL
        assocR
        subseg_cnt = []
        validity 
        max_err
        valid_flag = false
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
            obj.DataAssociation();% --> initial data association is done using init segmentation information
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
                obj.arc_segments{i}.L = ((obj.arc_segments{i}.L).^(1/2) + arc_delta_params{i}.L).^2;
                obj.arc_segments{i}.x0 = obj.arc_segments{i}.x0 + arc_delta_params{i}.x0;
                obj.arc_segments{i}.y0 = obj.arc_segments{i}.y0 + arc_delta_params{i}.y0;
                obj.arc_segments{i}.tau0 = obj.arc_segments{i}.tau0 + arc_delta_params{i}.tau0;
                obj.subseg_cnt(i) = length(arc_delta_params{i}.kappa);
            end

            % Perform data association with updated variables
%             obj.associate(); % update data association after one batch
%             optimization process for stable convergence
        end

        %% Map Visualization (2D Segmentwise Point Map)
        function obj = visualize2DMap(obj)
            % Alter this function for debugging purpose
            % Current : Visualize initial arc segmentation results

            figure(1); hold on; grid on; axis equal;
            
            n = length(obj.states);
            
            for i=1:n
                P = obj.states{i}.P;
                p_est = plot(P(1),P(2),'r.');
            end

            n = length(obj.arc_segments);
            for i=1:n
                m = length(obj.arc_segments{i}.kappa);
                seg = obj.arc_segments{i};
                heading = seg.tau0;
                SegPoints = [seg.x0;
                             seg.y0];
                for j=1:m
                    kappa = seg.kappa(j); L = seg.L(j);
                    headingPrev = heading;
                    heading = heading + kappa * L;
                    headingCurr = heading;

                    heading_ = linspace(headingPrev,headingCurr,1e3);
                    addedSegPoints = SegPoints(:,end) + 1/kappa * [sin(heading_) - sin(headingPrev);
                                                                   -cos(heading_) + cos(headingPrev)];
                    SegPoints = [SegPoints addedSegPoints];
                    
                    p_node = plot(addedSegPoints(1,1),addedSegPoints(2,1),'co');
                    plot(addedSegPoints(1,end),addedSegPoints(2,end),'co');
                end
                p_lane = plot(SegPoints(1,:),SegPoints(2,:),'k-');
                p_seg = plot(SegPoints(1,1),SegPoints(2,1),'ms');
                plot(SegPoints(1,end),SegPoints(2,end),'ms');

            end
           
            xlabel('Global X(m)'); ylabel('Global Y(m)');
            title('Optimized Vehicle Trajectory and Lane Segments');
            legend([p_est,p_lane,p_node,p_seg], ...
                   'Optimized Vehicle Trajectory', ...
                   'Optimized Arc Spline based ego-lane', ...
                   'Lane Sub-segment','Lane Segment');
        end
        
        %% Validate Current Map after convergence
        function obj = validate(obj,phase)
            % For each segment, find the most invalid subsegment
            n = size(obj.lane.FactorValidIntvs,1);
            
            % Validity: stores number of invalid measurements for each
            % sub-segment for every segment
            obj.validity = {};
            obj.max_err = {};
            for i=1:length(obj.arc_segments)
                obj.validity = [obj.validity {zeros(1,length(obj.arc_segments{i}.kappa))}];
                obj.max_err = [obj.max_err {zeros(1,length(obj.arc_segments{i}.kappa))}];
            end
            
            if phase == 2
                prev_max = 1;
            elseif phase == 3
                prev_max = obj.lane.prev_num;
            end

            for i=1:n
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                
                leftSegIdx = obj.segment_info(1,i);
                rightSegIdx = obj.segment_info(2,i);

                for j=lb:ub
                    for k=1:prev_max
                        leftSubSegIdx = obj.assocL(j,k);
                        rightSubSegIdx = obj.assocR(j,k);

                        if leftSubSegIdx > 0
                            [chisq,valid] = obj.isValid(leftSegIdx,leftSubSegIdx,j,k,'left');
                            if ~valid
                                obj.validity{leftSegIdx}(leftSubSegIdx) = obj.validity{leftSegIdx}(leftSubSegIdx) + 1;
                                obj.max_err{leftSegIdx}(leftSubSegIdx) = chisq;
                            end
                        end

                        if rightSubSegIdx > 0
                            [chisq,valid] = obj.isValid(rightSegIdx,rightSubSegIdx,j,k,'right');
                            if ~valid
                                obj.validity{rightSegIdx}(rightSubSegIdx) = obj.validity{rightSegIdx}(rightSubSegIdx) + 1;
                                obj.max_err{rightSegIdx}(rightSubSegIdx) = chisq;
                            end
                        end
                    end
                end
            end

            % Add new segments for the most "invalid" sub-segment  
            % if there are more than 2 "invalid" measurements
            
            disp('====================<Adding Segments>====================')
            obj.valid_flag = true;
            for i=1:length(obj.validity)
                % Choose Sub-segment with the most number of invalid
                % measurements as the replication target
                [max_val,badSubSegIdx] = max(obj.validity{i});
                if max_val > 2
                    obj.replicate(i,badSubSegIdx);
                    obj.valid_flag = false;
                end

%                 org_idxs = 1:length(obj.validity{i});
%                 trimmed_idxs = org_idxs(obj.validity{i} > 2);
%                 max_errs = obj.max_err{i}(obj.validity{i} > 2);
%                 
%                 if ~isempty(max_errs)                    
%                     [~,idx] = max(max_errs);
%                     badSubSegIdx = trimmed_idxs(idx);
%                     obj.replicate(i,badSubSegIdx);
%                     obj.valid_flag = false;
%                 end                
            end
            if obj.valid_flag
                disp('All segments are valid, optimization finished...')
            end
            disp('=========================================================')

            obj.DataAssociation();
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
                    
                % Sub-segment information
                for j=1:length(initSegments)
                    initParams.bnds = [initParams.bnds; initSegments{j}.bnds];
                    initParams.kappa = [initParams.kappa initSegments{j}.kappa];
                    initParams.L = [initParams.L initSegments{j}.L];
%                     if j == 1
%                         initParams.xc = initParams.x0 - 1/initParams.kappa(end) * sin(heading);
%                         initParams.yc = initParams.y0 + 1/initParams.kappa(end) * cos(heading);
%                     else
%                         initParams.xc = [initParams.xc initParams.xc(end) + (1/initParams.kappa(end-1) - 1/initParams.kappa(end)) * sin(heading)];
%                         initParams.yc = [initParams.yc initParams.yc(end) - (1/initParams.kappa(end-1) - 1/initParams.kappa(end)) * cos(heading)];
%                     end
                    heading = heading + initParams.kappa(end) * initParams.L(end);
                end
                % Find corresponding state_idx for initParams.bnds
                [~,c] = find(obj.segment_info == i);
                stateIntvs = obj.lane.FactorValidIntvs(c(1):c(end),1:2);
                initParams.stateIdxs = zeros(1,size(initParams.bnds,1));

                for j=1:size(initParams.bnds,1)
                    initParams.stateIdxs(j) = obj.findStateIdx(stateIntvs,initParams.bnds(j,2));
                end

                obj.arc_segments = [obj.arc_segments {initParams}];
                obj.subseg_cnt = [obj.subseg_cnt length(initParams.kappa)];                
            end
            

        end
        
        %% Initial Data Association
        function obj = DataAssociation(obj)
            % Initial Data Association for phase 2 optimization
            obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
            obj.assocR = obj.assocL;

            for i=1:size(obj.lane.FactorValidIntvs,1)
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);

                leftSegIdx = obj.segment_info(1,i);
                rightSegIdx = obj.segment_info(2,i);

                for j=lb:ub
                    % do not compare with rel idx
                    obj.assocL(j,1) = obj.DataAssocMatch(leftSegIdx,j);
                    obj.assocR(j,1) = obj.DataAssocMatch(rightSegIdx,j);
                end
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
            line_acc_thres = 0.1;

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
            rmse_thres = 0.5;
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
                
                err_thres = 1.5; 
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
        function obj = associate(obj,phase)
            % Data Association with current vehicle state and arc-spline
            % parameters. After each step in optimization, vehicle states
            % and arc parameters should be updated before re-association
            %
            % Create arc segment node points before data association
            % Data Association for phase 3 optimization
            n = length(obj.arc_segments);
            for i=1:n
                m = length(obj.arc_segments{i}.kappa);
                seg = obj.arc_segments{i};
                heading = seg.tau0;
                segPoints = [seg.x0;
                             seg.y0];
                for j=1:m            
                    % Propagate and save subsegment node points
                    kappa = seg.kappa(j); L = seg.L(j);
                    headingPrev = heading;
                    heading = heading + kappa * L;
                    headingCurr = heading;
                    segPoints = [segPoints segPoints(:,end) + 1/kappa * [sin(headingCurr) - sin(headingPrev);
                                                                        -cos(headingCurr) + cos(headingPrev)]];
                    
                    % Propagate and save subsegment arc center points
                    if j == 1
                        cPoints = [seg.x0 - 1/kappa * sin(headingPrev);
                                   seg.y0 + 1/kappa * cos(headingPrev)];
                    else
                        cPoints = [cPoints cPoints(:,end) + (1/seg.kappa(j-1) - 1/seg.kappa(j)) * [sin(headingPrev);-cos(headingPrev)]];
                    end
                end
                obj.arc_nodes = [obj.arc_nodes {segPoints}];
                obj.arc_centers = [obj.arc_centers {cPoints}];
            end

            % Data Association
            n = size(obj.lane.FactorValidIntvs,1);
            
            obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
            obj.assocR = obj.assocL;
            
            if phase == 2
                prev_max = 1;
            elseif phase == 3
                prev_max = obj.lane.prev_num;
            end
            
            for i=1:n
                % State lb and ub for each large segment
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);

                leftSegmentIdx = obj.segment_info(1,i);
                rightSegmentIdx = obj.segment_info(2,i);

                % Perform matching for valid points
                for j=lb:ub
                    lane_idx = obj.lane.state_idxs(j);                    
                    R = obj.states{j}.R;
                    P = obj.states{j}.P;
                    % Filter Valid state idx
                    if obj.lane.prob(lane_idx,2) > obj.lane_prob_thres && obj.lane.prob(lane_idx,3) > obj.lane_prob_thres
                        for k=1:prev_max
                            rel_pos = [10*(k-1);0;0];
                            pos = P + R * rel_pos;
                            
                            segIdxL = obj.match(pos,R,leftSegmentIdx);
                            segIdxR = obj.match(pos,R,rightSegmentIdx);

                            obj.assocL(j,k) = segIdxL;
                            obj.assocR(j,k) = segIdxR;
                        end
                    else % For low reliability lane measurement, use only 0m preview measurement                        
                        segIdxL = obj.match(P,R,leftSegmentIdx);
                        segIdxR = obj.match(P,R,rightSegmentIdx);

                        obj.assocL(j,1) = segIdxL;
                        obj.assocR(j,1) = segIdxR;
                    end
                end
            end
        end
        
        %% Check validity for each specific cases
        function [chisq,flag] = isValid(obj,segIdx,subsegIdx,state_idx,prev_idx,dir)
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
        function obj = replicate(obj,SegIdx,SubSegIdx)
            seg = obj.arc_segments{SegIdx};
            kappa = seg.kappa(SubSegIdx); L = seg.L(SubSegIdx);
            
            % We assume that for segment with one sub-segments are
            % well-optimized(Naively replicating this case will cause
            % singularity errors)

            if length(seg.kappa) == 1
                error(['Cannot replicate from segments with only one sub-segment,' ...
                       'lower initial parametrization error threshold to create' ...
                       'more than one sub-segment'])
            end

            % Create replica by halving the original arc length
            % Re-calculate curvature values
            % Update bnds values
            if SubSegIdx == 1
                seg.L(SubSegIdx) = 1/2 * seg.L(SubSegIdx);
                seg.L = [L/2 seg.L];
                
                rem_kappa = seg.kappa(SubSegIdx+1:end);
                next_kappa = seg.kappa(SubSegIdx+1);
                
%                 seg.kappa = [kappa, (kappa + next_kappa)/2, rem_kappa];   

%                 new_kappa = 2/(1/kappa + 1/next_kappa);
%                 seg.kappa = [kappa, new_kappa, rem_kappa];

                seg.kappa = [kappa, kappa, rem_kappa];
                
                rem_bnds = seg.bnds(SubSegIdx+1:end,:);
                curr_lb = seg.bnds(SubSegIdx,1); curr_ub = seg.bnds(SubSegIdx,2);
                new_bnd = ceil((curr_lb + curr_ub)/2);
                new_bnds = [curr_lb new_bnd;
                            new_bnd curr_ub];
                seg.bnds = [new_bnds; rem_bnds];

            elseif SubSegIdx == length(seg.kappa)                
                seg.L(SubSegIdx) = 1/2 * seg.L(SubSegIdx);
                seg.L = [seg.L L/2];
                rem_kappa = seg.kappa(1:SubSegIdx-1);
                prev_kappa = seg.kappa(SubSegIdx-1);
%                 seg.kappa = [rem_kappa, (prev_kappa + kappa)/2, kappa];
                
%                 new_kappa = 2/(1/kappa + 1/prev_kappa);
%                 seg.kappa = [rem_kappa, new_kappa, kappa];

                seg.kappa = [rem_kappa, kappa, kappa];

                rem_bnds = seg.bnds(1:SubSegIdx-1,:);
                curr_lb = seg.bnds(SubSegIdx,1); curr_ub = seg.bnds(SubSegIdx,2);
                new_bnd = ceil((curr_lb + curr_ub)/2);
                new_bnds = [curr_lb new_bnd;
                            new_bnd curr_ub];
                seg.bnds = [rem_bnds; new_bnds];

            else
                kappaF = seg.kappa(1:SubSegIdx-1);
                kappaB = seg.kappa(SubSegIdx+1:end);
                LF = seg.L(1:SubSegIdx-1);
                LB = seg.L(SubSegIdx+1:end);

%                 kappaF_ = (kappaF(end)+kappa)/2;
%                 kappaB_ = (kappa + kappaB(1))/2;
%                 seg.kappa = [kappaF kappaF_ kappaB_ kappaB];
                
%                 kappaF_ = 2/(1/kappaF(end) + 1/kappa);
%                 kappaB_ = 2/(1/kappa + 1/kappaB(1));
%                 seg.kappa = [kappaF kappaF_ kappaB_ kappaB];

                seg.kappa = [kappaF, kappa, kappa, kappaB];
                seg.L = [LF 1/2*L 1/2*L LB];

                bndsF = seg.bnds(1:SubSegIdx-1,:);
                bndsB = seg.bnds(SubSegIdx+1:end,:);
                curr_lb = seg.bnds(SubSegIdx,1); curr_ub = seg.bnds(SubSegIdx,2);
                new_bnd = ceil((curr_lb + curr_ub)/2);
                new_bnds = [curr_lb new_bnd;
                            new_bnd curr_ub];
                seg.bnds = [bndsF; new_bnds; bndsB];
            end
            % Update stateIdxs for segments
            [~,c] = find(obj.segment_info == SegIdx);
            stateIntvs = obj.lane.FactorValidIntvs(c(1):c(end),1:2);
            seg.stateIdxs = zeros(1,size(seg.bnds,1));
            for i=1:size(seg.bnds,1)
                seg.stateIdxs(i) = obj.findStateIdx(stateIntvs,seg.bnds(i,2));
            end

            % Update segment info
            obj.arc_segments{SegIdx} = seg;
            obj.subseg_cnt(SegIdx) = obj.subseg_cnt(SegIdx) + 1;
            disp(['SubSegment added at Segment Idx ',num2str(SegIdx),', SubSegment Idx ',num2str(SubSegIdx)])
        end
        
         %% Find initial data association match
         function idx = DataAssocMatch(obj,SegIdx,state_idx)
            idx = 0;
            [~,c] = find(obj.segment_info == SegIdx);
            stateIntv = obj.lane.FactorValidIntvs(c(1):c(end),:);
            diff = zeros(1,length(c));
            for i=1:length(c)-1
                diff(i+1) = stateIntv(i+1,1) - stateIntv(i,2) - 1;
            end
            summed_diff = cumsum(diff);
            
            rel_pos = 0;
            for i=1:length(c)
                if stateIntv(i,1) <= state_idx && stateIntv(i,2) >= state_idx
                    rel_pos = state_idx - stateIntv(1,1) + 1 - summed_diff(i);
                end
            end
            if rel_pos == 0
%                 disp(state_idx)
%                 disp(stateIntv)
%                 disp(diff)
                error('No appropriate match 1 ... Most likely to be implementation error')
            end

            Seg = obj.arc_segments{SegIdx};
            for i=1:size(Seg.bnds,1)
                if Seg.bnds(i,1) <= rel_pos && Seg.bnds(i,2) >= rel_pos
                    idx = i;
                end
            end
            if idx == 0
%                 disp(state_idx)
%                 disp(stateIntv)
%                 disp(diff)
%                 disp(rel_pos)
%                 disp(Seg.bnds)
                error('No appropriate match 2... Most likely to be implementation error')
            end
        end
        
        %% Point-SubSegment Matching
        function SubSegIdx = match(obj,pos,att,SegIdx)
            % May need to fix
            % Matching is performed in 2D manner
            segPoints = obj.arc_nodes{SegIdx};
            segCenters = obj.arc_centers{SegIdx};
            seg = obj.arc_segments{SegIdx};

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
                    dist = abs(sqrt(([xv;yv] - segCenters(:,i))' * ([xv;yv] - segCenters(:,i))) - abs(1/seg.kappa(i)));
                    
                    cand = [cand [i;dist]];
                end
            end

            if ~isempty(cand)
                % If match candidate exists, pick the one with smallest
                % matching error as the correct match
                [~,idx] = min(cand(2,:));
                SubSegIdx = cand(1,idx);
            else
                SubSegIdx = 0;
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
        
        %% Compute center points of arcs, given parameters
        function [xc,yc] = getCenter(x0,y0,tau0,kappa,L)
            heading = tau0;
            xc = []; yc = [];
            for j=1:length(kappa)
                if j == 1
                    xc = x0 - 1/kappa(j) * sin(heading);
                    yc = y0 + 1/kappa(j) * cos(heading);
                else
                    xc = xc + (1/kappa(j-1) - 1/kappa(j)) * sin(heading);
                    yc = yc - (1/kappa(j-1) - 1/kappa(j)) * cos(heading);
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
        
         %% Find state_idx from segment bnds
         function state_idx = findStateIdx(Intvs,num)
            augIdxs = [];            
            for i=1:size(Intvs,1)
                augIdxs = [augIdxs Intvs(i,1):Intvs(i,2)];
            end
            state_idx = augIdxs(num);
        end

    end
end