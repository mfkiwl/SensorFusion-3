classdef ArcMap < handle
%ARCMAP - Arc Spline based lane model for SNLS optimization
% 
%   When optimization mode is updated to '2-phase' in optimizer.m, this
%   module is automatically activated. Below is the overall framework 
%   of this module.
%   
%   Input: 
%   1. Optimized vehicle state (INS + GNSS + WSS Fusion)
%
%   2. Lane Measurement Information
%
%   ===================================================================
%   Step 1. Compute lane points using lane measurements and optimized
%           vehicle trajectory. Perform point clustering that should be
%           on the same lane.
%   
%   Step 2. Using divide & conquer method, find intervals for piecewise
%           linear approximation for point cluster obtained at step 1.
%           Then, merge linear approximations into arcs with some fixed
%           error threshold.
%           
%   Step 3. Using ArcFit.m, perform NLS optimization to find optimal
%           parameter values for the given number of arc parameters in
%           step 2.
%   
%   [From here, Mixed Usage with optimizer.m]
% 
%   Step 4. Use the results from step 3 as the initial arc spline map
%           and perform optimization using previewed lane measurements
%           in optimizer.m
% 
%   Step 5. Every iteration, perform data association(which lane 
%           measurement belongs to which Sub-segment). 
%
%   Step 6. Fully converge with fixed number of Sub-segments. Validate
%           current optimization results with lane measurement
%           reliability data. 
% 
%   Step 7. If invalid, add new segment to the current map and repeat
%           steps 4 ~ 6. If valid, end optimization.
%   ===================================================================
%   
%   Implemented by JinHwan Jeon, 2022

    properties (Access = public)
        states
        lane
        segments = {}
        LeftCov
        RightCov
        covs = {}
        ext_segments
        ext_covs = {}
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
        cubicFit = {}
        arcFit = {}
    end
    
    %% Public Methods
    methods (Access = public)
        %% Constructor
        function obj = ArcMap(states,lane,lane_prob_thres,LeftCov,RightCov)
            obj.states = states;
            obj.lane = lane;
            obj.lane_prob_thres = lane_prob_thres;
            obj.LeftCov = LeftCov;
            obj.RightCov = RightCov;
            obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
            obj.assocR = obj.assocL;

            % Initial Segment Classification
            obj.InitSegmentClassification();
            
%             % Initial Data Points Propagation
            obj.InitSegmentExtension();

            % Initial Segment Parametrization
            obj.InitSegmentParametrization();            
            

            % Perform Arc Spline based lane fitting with initial
            % parametrization data
%             for i=1:length(obj.arc_segments)
%                 % Fix from here! add covariance
%                 obj.dummy.initFit = ArcFit(obj.arc_segments{i},obj.segments{i},obj.arc_covs{i},obj.marker{i},i);    
%                 obj.dummy.initFit.optimize();
%                 
%                 obj.arc_segments{i} = obj.dummy.initFit.getParams();                
%             end
%             obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
%             obj.assocR = obj.assocL;
%             
%             obj.DataAssociation0();
%             obj.DataAssociationRem();     
            
            
            obj.cubicFit = {};
            % Perform Cubic Fitting (Comparison with Arc Splines)
            numSubSeg = 0;
            disp('[Cubic Spline Approximation Results]')
            disp('==============================================')
            for i=1:length(obj.segments)
                
                fit = CubicFit(obj.segments{i},obj.covs{i},i);
                fit.optimize();
                numSubSeg = numSubSeg + length(fit.segments);
                disp(['Segment ',num2str(i),': ',num2str(length(fit.segments)),' SubSegments'])
                obj.cubicFit = [obj.cubicFit, {fit}];
                fit.visualize();
            end
            disp('==============================================')
            disp(['Total number of Segments: ',num2str(length(obj.segments))])
            disp(['Total number of SubSegments: ',num2str(numSubSeg)])
            
            % Arc Fitting tested with EKF : Too many segments (Deprecated)
%             obj.arcFit = {};
%             % Perform Arc Fitting (Non Continuous)
%             numSubSeg = 0;
%             disp('  ')
%             disp('[Arc Spline Approximation Results]')
%             disp('==============================================')
%             for i=1:length(obj.segments)
%                 
%                 fit = ArcFitEKF(obj.segments{i},obj.covs{i},i);
%                 fit.optimize();
%                 numSubSeg = numSubSeg + length(fit.segments);
%                 disp(['Segment ',num2str(i),': ',num2str(length(fit.segments)),' SubSegments'])
%                 obj.arcFit = [obj.arcFit, {fit}];
%                 fit.visualize();
%             end
%             disp('==============================================')
%             disp(['Total number of Segments: ',num2str(length(obj.segments))])
%             disp(['Total number of SubSegments: ',num2str(numSubSeg)])
        end    
        
        %% Running Optimization for each Large Segment
        function obj = dummyF(obj)
            
            obj.arcFit = {};
            for i=4:length(obj.segments)
                
                initFit = ArcFit(obj.arc_segments{i}, ...
                                 obj.segments{i}, ...                                            
                                 obj.covs{i},i);    
                initFit.optimize();
                obj.arcFit = [obj.arcFit, {initFit}];
%                 obj.arc_segments{i} = obj.dummy.initFit.getParams();
                error('1')
            end
           
            
        end
        
        %% Optimize Parameters
        function obj = optimize(obj)

            for i=1:length(obj.arc_segments)
                initFit = ArcFit(obj.arc_segments{i}, ...
                                 obj.segments{i}, ...                                            
                                 obj.covs{i},i);                 
                obj.arc_segments{i} = initFit;
                obj.arc_segments{i}.optimize();
            end
%             obj.arcFit{1}.CreateTestArcFitBase();
        end

        %% Dummy func 3 for test
        function obj = dummy3(obj)
%            obj.arcFit{1}.optimize();
            obj.arcFit{1}.RunTestArcFitBase();
        end


        %% Map Update for 
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
            minL = obj.lane.minL;
            for i=1:length(obj.arc_segments)                
                obj.arc_segments{i}.kappa = obj.arc_segments{i}.kappa + arc_delta_params{i}.kappa;
                obj.arc_segments{i}.L = ((obj.arc_segments{i}.L - minL).^(1/2) + arc_delta_params{i}.L).^2 + minL;
                obj.arc_segments{i}.x0 = obj.arc_segments{i}.x0 + arc_delta_params{i}.x0;
                obj.arc_segments{i}.y0 = obj.arc_segments{i}.y0 + arc_delta_params{i}.y0;
                obj.arc_segments{i}.tau0 = obj.arc_segments{i}.tau0 + arc_delta_params{i}.tau0;
                obj.subseg_cnt(i) = length(arc_delta_params{i}.kappa);
            end

            % Perform data association with updated variables
            % Update data association for remaining preview distance 
            % measurements iteratively
            obj.assocL(:,2:obj.lane.prev_num) = zeros(size(obj.assocL,1),size(obj.assocL,2)-1);
            obj.assocR(:,2:obj.lane.prev_num) = zeros(size(obj.assocR,1),size(obj.assocR,2)-1);
            obj.DataAssociationRem();     
            % 0m previewed points have fixed data association 
            % changes only after validation(replicated segments)
        end
        
        %% Validate Current Map after convergence
        function obj = validate(obj)
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
            
            for i=1:n
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                
                leftSegIdx = obj.segment_info(1,i);
                rightSegIdx = obj.segment_info(2,i);

                for j=lb:ub
                    for k=1:obj.lane.prev_num
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
            end
            if obj.valid_flag
                % Finish optimization and display summary
                disp('All segments are valid, optimization finished...')
                obj.summary();
            end
            disp('=========================================================')
            
            % Perform Data Association for next optimization 
            obj.assocL = zeros(length(obj.states),obj.lane.prev_num);
            obj.assocR = obj.assocL;
            
            obj.DataAssociation0();
            obj.DataAssociationRem();
        end

        %% Map Visualization (2D Segmentwise Point Map)
        function obj = visualize2DMap(obj)
            % Visualize Parametrization Results            
                
            obj.summary();


            figure(50); hold on; grid on; axis equal;
            
            n = length(obj.states);
            
            for i=1:n
                P = obj.states{i}.P;
                p_est = plot(P(1),P(2),'r.');
            end

            n = length(obj.arc_segments);
            % 1:n
            for i=1:n
                m = length(obj.arc_segments{i}.params.kappa);
                seg = obj.arc_segments{i}.params;
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
        
        %% Display Optimization Summary
        function summary(obj)
            disp('[Summary of Optimized Arc Segments and Subsegments]')
            
            numSubSeg = 0;
            for i=1:length(obj.arc_segments)
                seg = obj.arc_segments{i}.params;
                numSubSeg = numSubSeg + length(seg.L);
                fprintf('\n')
                disp(['Segment No. ',num2str(i),': ',num2str(length(seg.L)),' Subsegments'])
                
                for j=1:length(seg.L)
                    disp(['   Subsegment No. ',num2str(j), ...
                          ': L = ',num2str(seg.L(j),3), ...
                          'm, R = ',num2str(abs(1/seg.kappa(j)),3),'m'])
                end
                fprintf('\n')
                disp(['   Segment Length: ',num2str(sum(seg.L),3),'m'])
            end
            fprintf('\n')
            disp('[Summary]')
            disp(['Number of Segments: ',num2str(length(obj.arc_segments)), ...
                  ', Number of Subsegments: ',num2str(numSubSeg)])
        end   
        
        %% Check if data association has gone wrong
        function run(obj)
            obj.InitSegmentExtension();
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
        end
        
        %% Initial Segment Extension
        function obj = InitSegmentExtension(obj)
            % Propagate lane points based on raw lane measurements
            % segments: 0m previewed lane points
            % ext_segments: remaining previewed lane points
            % 

            obj.segments = {};
            obj.covs = {};
%             obj.ext_segments = {}; % initialize segments variable again
%             obj.ext_covs = {};
            n = size(obj.lane.FactorValidIntvs,1);            
            
            for i=1:n+1
                obj.segments = [obj.segments {[]}];
                obj.covs = [obj.covs {[]}];
%                 obj.ext_segments = [obj.ext_segments {[]}];
%                 obj.ext_covs = [obj.ext_covs {[]}];                
            end

            for i=1:n
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                Left0 = []; Right0 = []; LeftCov0 = []; RightCov0 = [];
%                 Left = []; Right = []; LeftCovR = []; RightCovR = [];                
            
                for j=lb:ub
                    R = obj.states{j}.R;
                    P = obj.states{j}.P;
%                     rpy = dcm2rpy(R); psi = rpy(3);
%                     R2d = [cos(psi), -sin(psi);
%                            sin(psi), cos(psi)];
                    
%                     lane_idx = obj.lane.state_idxs(j);
%                     Lprob = obj.lane.prob(lane_idx,2);
%                     Rprob = obj.lane.prob(lane_idx,3);
                    
                    LeftL = P + R * obj.states{j}.left(:,1);
                    RightL = P + R * obj.states{j}.right(:,1);
                    Left0 = [Left0, LeftL(1:2)];
                    Right0 = [Right0, RightL(1:2)];
                    LeftCov0 = [LeftCov0, obj.LeftCov(4*j-3:4*j,1)];
                    RightCov0 = [RightCov0, obj.RightCov(4*j-3:4*j,1)];                    
                    
%                     % Add more lane measurements 
%                     if Lprob > obj.lane_prob_thres && Rprob > obj.lane_prob_thres
% 
%                         for k=2:obj.lane.prev_num
%                             Lstd = obj.lane.lystd(lane_idx,k);
%                             Rstd = obj.lane.rystd(lane_idx,k);
%                             
%                             LeftL = P + R * obj.states{j}.left(:,k);
%                             RightL = P + R * obj.states{j}.right(:,k);
%                             if Lstd < 3.4 && Rstd < 3.4
%                                 Left = [Left, LeftL(1:2)];
%                                 Right = [Right, RightL(1:2)];
%                                 LeftCovR = [LeftCovR, obj.LeftCov(4*j-3:4*j,k)];
%                                 RightCovR = [RightCovR, obj.RightCov(4*j-3:4*j,k)];
%                             end
%                         end                        
%                     end                    
                end                
                
                obj.segments{obj.segment_info(1,i)} = [obj.segments{obj.segment_info(1,i)} Left0];
                obj.segments{obj.segment_info(2,i)} = [obj.segments{obj.segment_info(2,i)} Right0];
                obj.covs{obj.segment_info(1,i)} = [obj.covs{obj.segment_info(1,i)} LeftCov0];
                obj.covs{obj.segment_info(2,i)} = [obj.covs{obj.segment_info(2,i)} RightCov0];

%                 obj.ext_segments{obj.segment_info(1,i)} = [obj.ext_segments{obj.segment_info(1,i)} Left];
%                 obj.ext_segments{obj.segment_info(2,i)} = [obj.ext_segments{obj.segment_info(2,i)} Right];
%                 obj.ext_covs{obj.segment_info(1,i)} = [obj.ext_covs{obj.segment_info(1,i)} LeftCovR];
%                 obj.ext_covs{obj.segment_info(2,i)} = [obj.ext_covs{obj.segment_info(2,i)} RightCovR];
            end

%             disp('[Eliminating invalid extensions...]')
%             % Eliminate points that are outside of the region of interest
%             for i=1:length(obj.segments) 
%                 seg = obj.segments{i};
%                 ext_seg = obj.ext_segments{i};
%                 org_n = size(ext_seg,2);
%                 bools = [];
%                 % 
%                 valid = false;
%                 j = 0;
% %                 if i == 10
% %                     figure(1);hold on; grid on; axis equal;
% %                     plot(seg(1,:),seg(2,:),'r.'); 
% %                     plot(ext_seg(1,:),ext_seg(2,:),'b.');
% %                     disp(length(ext_seg))
% %                 end
%                 if org_n > 0
%                     while ~valid
%                         j = j + 1;
% 
%                         % May need to change 10 to other larger values for
%                         % more accurate removal
%                         valid = obj.isMatched(seg(:,end-10:end),ext_seg(:,j)); 
%                         if valid
%                             flag = true;
%                             break;
%                         else
%                             bools = [bools true];
%                             
%                             if j == org_n % No Valid Match (Too short segment)
%                                 obj.ext_segments{i} = [];
%                                 new_n = size(obj.ext_segments{i},2);
%                                 disp(['[Segment ',num2str(i),'] : ',num2str(org_n-new_n),' points discarded'])
%                                 flag = false;
%                                 break;
%                             end
%                         end                    
%                     end                
%                     
%                     if flag
%     %                     disp(num2str(i))
%                         while j <= length(ext_seg)
%                             indc = obj.isMatched(seg,ext_seg(:,j));
%                             bools = [bools indc];    
%                             j = j + 1;
%                         end
%                         
%                         obj.ext_segments{i} = ext_seg(:,bools == 1);
%                         new_n = size(obj.ext_segments{i},2);
%                         disp(['[Segment ',num2str(i),'] : ',num2str(org_n-new_n),' points discarded'])
%                     end
%                 else
%                     disp(['[Segment ',num2str(i),'] : Already Empty!'])
%                 end
%             end
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
            rmse_thres = 0.2;
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
                iter_cnt = 1;
                iter_lim = n - lb + 1;
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
                    iter_cnt = iter_cnt + 1;
                    % If Looped too many times, exit with ub such that
                    % error is larger than 0.5
                    if iter_cnt >= iter_lim
                        err_thres = 0.5; % Check if this is valid
                    end
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
        
        %% Data Association for 0m preview
        function obj = DataAssociation0(obj)
            n = length(obj.arc_segments);

            for i=1:n
                seg = obj.arc_segments{i};
                m = size(seg.bnds,1);
                [r,c] = find(obj.segment_info == i);
                stateIntvs = obj.lane.FactorValidIntvs(c(1):c(end),1:2);     

                augIdxs = [];      
                augDirs = [];
                for j=1:size(stateIntvs,1)                
                    augIdxs = [augIdxs stateIntvs(j,1):stateIntvs(j,2)];
                    augDirs = [augDirs r(j) * ones(1,stateIntvs(j,2)-stateIntvs(j,1)+1)];
                end

                for j=1:m
                    stateIdxs = augIdxs(seg.bnds(j,1):seg.bnds(j,2));
                    dirs = augDirs(seg.bnds(j,1):seg.bnds(j,2));                    
                    
                    l = length(stateIdxs);
                    for k=1:l
                        if dirs(k) == 1
                            obj.assocL(stateIdxs(k),1) = j;
                        elseif dirs(k) == 2
                            obj.assocR(stateIdxs(k),1) = j;
                        end
                    end
                end
            end
        end

        %% Data Association for remaining preview
        function obj = DataAssociationRem(obj)
            % Data Association with current vehicle state and arc-spline
            % parameters. After each step in optimization, vehicle states
            % and arc parameters should be updated before re-association
            %
            % Pre-compute arc segment node/center points 
            
            n = length(obj.arc_segments);
            obj.arc_nodes = {};
            obj.arc_centers = {};
            for i=1:n                
                seg = obj.arc_segments{i};
                initParams = [seg.x0, seg.y0, seg.tau0];
                kappa = seg.kappa; L = seg.L;                
                obj.arc_nodes = [obj.arc_nodes {obj.propNode(initParams,kappa,L)}];
                obj.arc_centers = [obj.arc_centers {obj.propCenter(initParams,kappa,L)}];
            end

            % Data Association
            n = size(obj.lane.FactorValidIntvs,1);                        

            for i=1:n
                % State lb and ub for each large segment
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);

                leftSegIdx = obj.segment_info(1,i);
                rightSegIdx = obj.segment_info(2,i);

                % Perform matching for valid points
                for j=lb:ub
                    lane_idx = obj.lane.state_idxs(j);                    
                    R = obj.states{j}.R;
                    P = obj.states{j}.P;
                    
                    rpy = dcm2rpy(R); psi = rpy(3);
                    R2d = [cos(psi), -sin(psi);sin(psi), cos(psi)];

                    % Filter Valid state idx
                    if obj.lane.prob(lane_idx,2) > obj.lane_prob_thres && obj.lane.prob(lane_idx,3) > obj.lane_prob_thres
                        % Use preview distance measurement only for
                        % reliable measurements
                        
                        for k=2:obj.lane.prev_num    
                            % If std value of previewed measurement is too
                            % large(abnormal values), reject measurement
                            % even if the lane probability is measured to
                            % be valid. 
                            % 
                            % This is to prevent impossible lane
                            % measurements from corrupting optimization
                            if obj.lane.lystd(lane_idx,k) < 3.4 || obj.lane.rystd(lane_idx,k) < 3.4                                

                                rel_pos = [10*(k-1);0];
                                pos = P(1:2) + R2d * rel_pos;  
                                
                                segIdxL = obj.match(pos,R,leftSegIdx);
                                segIdxR = obj.match(pos,R,rightSegIdx);
                            
                                obj.assocL(j,k) = segIdxL;
                                obj.assocR(j,k) = segIdxR;                                
                            end                            
                        end
%                     else
%                         % Perform data association only for 0m preview
%                         % measurement                        
%                         segIdxL = obj.match(P,R,leftSegIdx);
%                         segIdxR = obj.match(P,R,rightSegIdx);
%                         obj.assocL(j,1) = segIdxL;
%                         obj.assocR(j,1) = segIdxR;
                    end

                    % Check for any abnormal zeros                    
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
            % * Update bnds values *
            if SubSegIdx == 1
                seg.L(SubSegIdx) = 1/2 * seg.L(SubSegIdx);
                seg.L = [L/2 seg.L];
%                 next_L = seg.L(SubSegIdx+1);
%                 curr_L = seg.L(SubSegIdx);
%                 new_L = 1/2 * (1/2 * curr_L + next_L);
%                 rem_L = seg.L(SubSegIdx+2:end);
%                 seg.L = [curr_L/2, new_L, new_L, rem_L];
                
                rem_kappa = seg.kappa(SubSegIdx+1:end);
%                 next_kappa = seg.kappa(SubSegIdx+1);
%                 new_kappa = 2/(1/kappa + 1/next_kappa);
%                 seg.kappa = [kappa, new_kappa, rem_kappa];

                seg.kappa = [kappa, kappa, rem_kappa];
                
                curr_bnd = seg.bnds(SegIdx,:);
                next_bnd = seg.bnds(SegIdx+1,:);
                
                div1 = floor((curr_bnd(1)+curr_bnd(2))/2);
                div2 = floor((div1 + next_bnd(2))/2);

                rem_bnds = seg.bnds(SegIdx+2:end,:);

                seg.bnds = [curr_bnd(1), div1;
                            div1, div2;
                            div2, next_bnd(2);
                            rem_bnds];

            elseif SubSegIdx == length(seg.kappa)                
                seg.L(SubSegIdx) = 1/2 * seg.L(SubSegIdx);
                seg.L = [seg.L L/2];
                
%                 prev_L = seg.L(SubSegIdx-1);
%                 curr_L = seg.L(SubSegIdx);
%                 new_L = 1/2 * (1/2 * curr_L + prev_L);
%                 rem_L = seg.L(1:SubSegIdx-2);                
%                 seg.L = [rem_L, new_L, new_L, curr_L/2];

                rem_kappa = seg.kappa(1:SubSegIdx-1);
%                 prev_kappa = seg.kappa(SubSegIdx-1);                
%                 new_kappa = 2/(1/kappa + 1/prev_kappa);
%                 seg.kappa = [rem_kappa, new_kappa, kappa];

                seg.kappa = [rem_kappa, kappa, kappa];

                curr_bnd = seg.bnds(SubSegIdx,:);
                prev_bnd = seg.bnds(SubSegIdx-1,:);

                div1 = floor((curr_bnd(1) + curr_bnd(2))/2);
                div2 = floor((prev_bnd(1) + div1)/2);

                rem_bnds = seg.bnds(1:SubSegIdx-2,:);
                
                seg.bnds = [rem_bnds;
                            prev_bnd(1), div2;
                            div2, div1;
                            div1, curr_bnd(2)];
            else
                kappaF = seg.kappa(1:SubSegIdx-1);
                kappaB = seg.kappa(SubSegIdx+1:end);
                LF = seg.L(1:SubSegIdx-1);
                LB = seg.L(SubSegIdx+1:end);
                curr_L = seg.L(SubSegIdx);
                seg.L = [LF, curr_L/2, curr_L/2, LB];

%                 next_L = seg.L(SubSegIdx+1);
%                 curr_L = seg.L(SubSegIdx);
%                 prev_L = seg.L(SubSegIdx-1);
% 
%                 new_LF = 1/2 * (prev_L + 1/2 * curr_L);                
%                 new_LB = 1/2 * (next_L + 1/2 * curr_L);
%                 seg.L = [LF, new_LF, new_LF, new_LB, new_LB, LB];
                

%                 kappaF_ = 2/(1/kappaF(end) + 1/kappa);
%                 kappaB_ = 2/(1/kappa + 1/kappaB(1));
%                 seg.kappa = [kappaF kappaF_ kappaB_ kappaB];

                seg.kappa = [kappaF, kappa, kappa, kappaB];

                curr_bnd = seg.bnds(SubSegIdx,:);
                prev_bnd = seg.bnds(SubSegIdx-1,:);
                next_bnd = seg.bnds(SubSegIdx+1,:);

                divC = floor((curr_bnd(1) + curr_bnd(2))/2);
                divP = floor((divC + prev_bnd(1))/2);
                divN = floor((divC + next_bnd(2))/2);

                P_bnds = seg.bnds(1:SubSegIdx-2,:);
                N_bnds = seg.bnds(SubSegIdx+2:end,:);

                seg.bnds = [P_bnds;
                            prev_bnd(1), divP;
                            divP, divC;
                            divC, divN;
                            divN, next_bnd(2);
                            N_bnds];
            end

            % Update segment info
            obj.arc_segments{SegIdx} = seg;
            obj.subseg_cnt(SegIdx) = obj.subseg_cnt(SegIdx) + 1;
            disp(['SubSegment added at Segment Idx ',num2str(SegIdx),', SubSegment Idx ',num2str(SubSegIdx)])
        end
        
        %% Point-SubSegment Matching for remaining preview measurement
        function SubSegIdx = match(obj,pos,att,SegIdx)
            % May need to fix
            % Matching is performed in 2D manner
            segNodes = obj.arc_nodes{SegIdx};
            
%             segCenters = obj.arc_centers{SegIdx};
%             seg = obj.arc_segments{SegIdx};

            n = size(segNodes,2);
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
                
                P1 = segNodes(:,i); P2 = segNodes(:,i+1);
                dP = (P1 - P2)' * (P1 - P2);
                
                x1 = P1(1); y1 = P1(2); 
                x2 = P2(1); y2 = P2(2);
                slope = (y2 - y1)/(x2 - x1);

                x_intersect = (x1 * (y2 - y1)/(x2 - x1) - y1 + xv * 1/tan(psi) + yv)/((y2 - y1)/(x2 - x1) + 1/tan(psi));
                y_intersect = slope * (x_intersect - x1) + y1; 
                P_intersect = [x_intersect; y_intersect];

                d1 = (P1 - P_intersect)' * (P1 - P_intersect);
                d2 = (P2 - P_intersect)' * (P2 - P_intersect);
                [~,idx] = max([dP,d1,d2]);
                if idx == 1
                    % Matching Success
                    % Compute lane point to intersection point distance
                    dist = (xv - x_intersect)^2 + (yv - y_intersect)^2;                      
                    cand = [cand [i;dist]];
                end
            end

            if ~isempty(cand)
                % If multiple match candidate exists, pick the one with
                % smallest matching error as the correct match
                % Vehicle to intersection point distance is small if they
                % are correctly matched
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
        
        %% Find Maximum column index
        function [idx, dir] = maxColIdx(arr,num)
            [r,c] = find(arr==num);
            if isempty(c)
                error('No matched number')
            else
                [idx,max_idx] = max(c);
                if r(max_idx) == 1
                    dir = 'left';
                elseif r(max_idx) == 2
                    dir = 'right';
                end
            end
        end

        %% Find Minimum column index
        function [idx, dir] = minColIdx(arr,num)
            [r,c] = find(arr==num);
            if isempty(c)
                error('No matched number')
            else
                [idx,min_idx] = min(c);
                if r(min_idx) == 1
                    dir = 'left';
                elseif r(min_idx) == 2
                    dir = 'right';
                end
            end
        end

        %% Find state_idx from segment bnds
        function [state_idxs,dirs] = findStateIdx(Intvs,num,dir_arr)
            augIdxs = [];      
            augDirs = [];
            for i=1:size(Intvs,1)                
                augIdxs = [augIdxs Intvs(i,1):Intvs(i,2)];
                augDirs = [augDirs dir_arr(i) * ones(1,Intvs(i,2)-Intvs(i,1)+1)];
            end
            state_idxs = augIdxs(num);
            dirs = augDirs(num);
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
        
        %% Reorder Lane Point Measurements
        function [reorderedLP, reorderedCov] = reorder(LP,Cov,endPoint)
            n = size(LP,2);
            LP = vertcat(LP,1:n);
            cpydLP = vertcat(LP(:,2:end),2:1:n); 
            cpydCov = Cov(:,2:end);
            

            reorderedLP = LP(1:2,1);
            reorderedCov = Cov(:,1);
            

            % Perform Nearest Neighbor Search
            search_idx = 1;
            d = sqrt((LP(1:2,search_idx) - endPoint)' * (LP(1:2,search_idx) - endPoint));
            while d > 0.5
                px = LP(1,search_idx); py = LP(2,search_idx);
                D = (cpydLP(1,:) - px).^2 + (cpydLP(2,:) - py).^2;

                [~,idx] = min(D);
                search_idx = cpydLP(3,idx);
                cpydLP(:,idx) = []; cpydCov(:,idx) = []; 
                reorderedLP = [reorderedLP LP(1:2,search_idx)];
                reorderedCov = [reorderedCov Cov(:,search_idx)];
                

                d = sqrt((LP(1:2,search_idx) - endPoint)' * (LP(1:2,search_idx) - endPoint));
            end

            if d > 1e-4
                reorderedLP = [reorderedLP endPoint];
                px = LP(1,:); py = LP(2,:);
                D = sqrt((px - endPoint(1)).^2 + (py - endPoint(2)).^2);
                [~,idx] = min(D);
                reorderedCov = [reorderedCov Cov(:,idx)];
                
            end
        end
        
        %% Return matching boolean
        function indc = isMatched(LP,Point)
            n = size(LP,2);
            indc = false;
            for i=1:n-1
                P1 = LP(:,i); P2 = LP(:,i+1);
                x1 = P1(1); y1 = P1(2); x2 = P2(1); y2 = P2(2);
                xp = Point(1); yp = Point(2);
                slope = (y2-y1)/(x2-x1);

                x_vert = (xp + x1 * slope^2 + (yp-y1)*slope)/(1 + slope^2);
                y_vert = slope * (x_vert - x1) + y1;
                X_vert = [x_vert; y_vert];

                Dp = (P1 - P2)' * (P1 - P2);
                D1 = (P1 - X_vert)' * (P1 - X_vert);
                D2 = (P2 - X_vert)' * (P2 - X_vert);
%                 disp(['Idx ',num2str(i),'~',num2str(i+1)])
%                 disp([Dp,D1,D2])
                [~,max_idx] = max([Dp,D1,D2]);  
                if max_idx == 1
                    dist = sqrt((xp - x_vert)^2 + (yp - y_vert)^2);
                    % If matched distance is too large, not a valid match
                    % This is to consider high-curvature road intervals,
                    % where multiple matches are possible
                    if dist < 5
                        indc = true;
                        break;
                    end
                end

            end
        end

    end
end