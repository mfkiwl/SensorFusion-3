classdef Map < handle
    % MAP - map class for lane factor based optimization
    % Initialized by creating virtual map of 0m previewed points
    % * Lane Changes are considered for creating segments
    % * For every segment, last states' previewed points are used to extend
    % virtual map
    % 
    % Optimization Scheme
    % 1. Perform association of lane points 
    % --> Optimize one step with that information
    % 2. Update association using retracted states 
    % 3. Repeat until convergence
    %
    % Implemented by JinHwan Jeon, 2022
    
    properties (Access = public)
        states
        lane
        segments = {}
        extended_segments = {}
        lane_prob_thres
        assoc = {}
    end

    %% Public Methods
    methods (Access = public)
        %% Constructor
        function obj = Map(states,lane,lane_prob_thres)
            obj.states = states;
            obj.lane = lane;
            obj.lane_prob_thres = lane_prob_thres;

            % Initialize Map 
            obj.initialize();

            % Association
            obj.associate();

        end

        %% Map Update
        function obj = update(obj)

        end

    end

    %% Private Methods
    methods (Access = private)
        %% Initialization
        function obj = initialize(obj)
            
            n = size(obj.lane.FactorValidIntvs,1);
            for i=1:n
                map = struct();
                map.left = [];
                map.right = [];
                map.will_be_matched = [];

                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2);
                
                % Except last state of the segment, add 0m previewed points
                for j = lb:ub
                    R = obj.states{j}.R;
                    P = obj.states{j}.P;
                    map.left = [map.left P + R * obj.states{j}.left(:,1)];
                    map.right = [map.right P + R * obj.states{j}.right(:,1)];
                    state_idx = obj.lane.state_idxs(j);

                    % If any of left/right lane detection probability is
                    % lower than threshold, no matching for that state
                    % index
                    if obj.lane.prob(state_idx,2) < obj.lane_prob_thres || obj.lane.prob(state_idx,3) < obj.lane_prob_thres
                        map.will_be_matched = [map.will_be_matched false];
                    else
                        map.will_be_matched = [map.will_be_matched true];
                    end
                end
                obj.segments = [obj.segments {map}];

                % Create Extended Map using last timestep previewed points
                
                R = obj.states{ub}.R;
                P = obj.states{ub}.P;
                ext_map = struct();
                ext_map.left = P + R * obj.states{ub}.left(:,2:end);
                ext_map.right = P + R * obj.states{ub}.right(:,2:end);
                
                obj.extended_segments = [obj.extended_segments {ext_map}];
            end
        end
        
        %% Lane Point data association 
        % Currently Fixing errors
        % 
        % Issue 1: Wrong match with interpolated measurement and actual
        % target match point. 
        % 
        %

        function obj = associate(obj)
            n = size(obj.lane.FactorValidIntvs,1);
            sampleXthres = 10; % may need to change threshold for valid interval sampling
            for i=1:n
                map = obj.MergeMap(obj.segments{i},obj.extended_segments{i});
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2); 
                assocL = zeros(ub - lb,obj.lane.prev_num-1);
                assocLprev = ones(ub - lb,obj.lane.prev_num-1);
                assocR = assocL;
                assocRprev = assocLprev;
                
                % Association for "normal" cases
                for j=lb:ub-1
                    if map.will_be_matched(j)
                        % Perform map matching only for states with good lane
                        % detection reliability
                        R = obj.states{j}.R;
                        P = obj.states{j}.P;
                        LEFT = obj.states{j}.left;
                        RIGHT = obj.states{j}.right;
                        % Perform frame transformation w.r.t vehicle frame
                        map_b = obj.frameTransformation(R,P,map);  
                        xL = map_b.left(1,:);
                        xR = map_b.right(1,:);

                        for k=2:obj.lane.prev_num
                            idxL = find(xL > 10*(k-1) - sampleXthres & xL < 10*(k-1) + sampleXthres); 
                            idxR = find(xR > 10*(k-1) - sampleXthres & xR < 10*(k-1) + sampleXthres);
                            disp(['Target State Idx: ', num2str(j)])
                            disp(['Preview Number: ', num2str(k)])
                            if isempty(idxL) || isempty(idxR)
                                error('No valid idx after transformation...')
                            end
                            candL = {}; candR = {};
                            % Left
                            for t=1:length(idxL)-1
                                idxP = idxL(t); idxN = idxL(t+1);
                                if xL(idxP) < 10*(k-1) && xL(idxN) >= 10*(k-1)
                                    % If match candidate is found,
                                    % interpolate and find resample point.
                                    % There may be multiple candidates
                                    
                                    disp([xL(idxP), xL(idxN)])

                                    [state_idxP,preview_idxP] = obj.getStateIdx(idxP,i);
                                    [state_idxN,preview_idxN] = obj.getStateIdx(idxN,i);

                                    disp([state_idxP,preview_idxP])
                                    disp([state_idxN,preview_idxN])

                                    Rt = obj.states{state_idxP}.R;
                                    Pt = obj.states{state_idxP}.P;

                                    lane_idxP = obj.lane.state_idxs(state_idxP);
                                    lane_idxN = obj.lane.state_idxs(state_idxN);
                                    x_ = [xL(idxP) xL(idxN)];
                                    stdy_ = [obj.lane.lystd(lane_idxP,preview_idxP),...
                                             obj.lane.lystd(lane_idxN,preview_idxN)];
                                    stdz_ = [obj.lane.lzstd(lane_idxP,preview_idxP),...
                                             obj.lane.lystd(lane_idxN,preview_idxN)];
                                    y_ = [obj.lane.ly(lane_idxP,preview_idxP),...
                                          obj.lane.ly(lane_idxN,preview_idxN)];
                                    z_ = [obj.lane.lz(lane_idxP,preview_idxP),...
                                          obj.lane.lz(lane_idxN,preview_idxN)];

                                    stdy_inter = interp1(x_,stdy_,10*(k-1));
                                    stdz_inter = interp1(x_,stdz_,10*(k-1));
                                    y_inter = interp1(x_,y_,10*(k-1));
                                    z_inter = interp1(x_,z_,10*(k-1));
                                    meas_inter = [y_inter;z_inter];
                                    cov_inter = diag([stdy_inter^2,stdz_inter^2]);

                                    cand = struct();
                                    [cand.idxP, cand.prev_idx] = obj.getStateIdx(idxP,i);
                                    cand.meas = meas_inter;
                                    cand.cov = cov_inter;
                                    % Need to convert actual world frame
                                    % lane point coords to matching frame

                                    cand.actual_meas = Rt' * (P + R * LEFT(:,k) - Pt);
                                    candL = [candL {cand}];
                                                    
                                    disp(['Interpolated Measurement: ', num2str(cand.meas')])
                                    disp(['Target Match Point: ', num2str(cand.actual_meas')])
                                end                                
                            end

                            % Right

                            for t=1:length(idxR)-1
                                idxP = idxR(t); idxN = idxR(t+1);
                                if xR(idxP) < 10*(k-1) && xR(idxN) >= 10*(k-1)
                                    % If match candidate is found,
                                    % interpolate and find resample point.
                                    % There may be multiple candidates
                                    
                                    [state_idxP,preview_idxP] = obj.getStateIdx(idxP,i);
                                    [state_idxN,preview_idxN] = obj.getStateIdx(idxN,i);
                                    Rt = obj.states{state_idxP}.R;
                                    Pt = obj.states{state_idxP}.P;

                                    lane_idxP = obj.lane.state_idxs(state_idxP);
                                    lane_idxN = obj.lane.state_idxs(state_idxN);
                                    x_ = [xR(idxP) xR(idxN)];
                                    stdy_ = [obj.lane.rystd(lane_idxP,preview_idxP),...
                                             obj.lane.rystd(lane_idxN,preview_idxN)];
                                    stdz_ = [obj.lane.rzstd(lane_idxP,preview_idxP),...
                                             obj.lane.rystd(lane_idxN,preview_idxN)];
                                    y_ = [obj.lane.ry(lane_idxP,preview_idxP),...
                                          obj.lane.ry(lane_idxN,preview_idxN)];
                                    z_ = [obj.lane.rz(lane_idxP,preview_idxP),...
                                          obj.lane.rz(lane_idxN,preview_idxN)];

                                    stdy_inter = interp1(x_,stdy_,10*(k-1));
                                    stdz_inter = interp1(x_,stdz_,10*(k-1));
                                    y_inter = interp1(x_,y_,10*(k-1));
                                    z_inter = interp1(x_,z_,10*(k-1));
                                    meas_inter = [y_inter;z_inter];
                                    cov_inter = diag([stdy_inter^2,stdz_inter^2]);

                                    cand = struct();
                                    [cand.idxP, cand.prev_idx] = obj.getStateIdx(idxP,i);
                                    cand.meas = meas_inter;
                                    cand.cov = cov_inter;
                                    % Need to convert actual world frame
                                    % lane point coords to matching frame

                                    cand.actual_meas = Rt' * (P + R * RIGHT(:,k) - Pt);
                                    candR = [candR {cand}];
                                end                                
                            end

                            % Verify candidates by performing Chi-Square
                            % Test. 
                            [validityL,NL] = obj.ChiSquareTest(candL);
                            [validityR,NR] = obj.ChiSquareTest(candR);

                            % Find acceptable candidate from validity array
                            % output valid index is the index of best
                            % candidate
                            validIdxL = obj.ProcValidity(validityL,NL);
                            validIdxR = obj.ProcValidity(validityR,NR);
                            
                            bestIdxL = candL{validIdxL}.idxP;
                            bestIdxL_prev = candL{validIdxL}.prev_idx;

                            bestIdxR = candR{validIdxR}.idxP;
                            bestIdxR_prev = candR{validIdxR}.prev_idx;

                            assocL(j-lb+1,k-1) = bestIdxL;
                            assocLprev(j-lb+1,k-1) = bestIdxL_prev;
                            assocR(j-lb+1,k-1) = bestIdxR;
                            assocRprev(j-lb+1,k-1) = bestIdxR_prev;
                        end

                        error('1')
                    end
                end
                
                % Save Association Matrix for one segment
                assoc_ = struct();
                assoc_.L = assocL;
                assoc_.R = assocR;
                assoc_.Lprev = assocLprev;
                assoc_.Rprev = assocRprev;
                obj.assoc = [obj.assoc {assoc_}];

                % Data Association for Extended Map Matching
                if i < n
                    
                
                end
            end

        end
        
        %% Get State number from map index and segment number
        function [state_idx,preview_idx] = getStateIdx(obj,map_idx,seg_idx)
            curr_seg_intv = obj.lane.FactorValidIntvs(seg_idx,1:2);
            bds = curr_seg_intv(1):curr_seg_intv(2);
            
            if map_idx > length(bds)
                % Current matched point is part of the extended map
                state_idx = bds(end);
                preview_idx = map_idx - length(bds) + 1;
            else
                % Current matched point is part of the normal map
                state_idx = bds(map_idx);
                preview_idx = 1;
            end
            
        end

    end

    %% Static Methods
    methods (Static)
        %% Convert map points into vehicle frame
        function map_b = frameTransformation(R,p,map)
            map_b = struct();
            map_b.left = R' * map.left - p;
            map_b.right = R' * map.right - p;
        end
        
        %% Perform Chi-Square Test for matching candidates
        function [validity,N] = ChiSquareTest(cand)
            n = length(cand);
            validity = false(1,n);
            N = zeros(1,n);
            % 99% Reliability for 2 DOF system
            thres = chi2inv(0.99,2); 

            for i=1:n
                act_meas = cand{i}.actual_meas(2:3);                
                meas = cand{i}.meas;
                cov = cand{i}.cov;
                
                % Compute Mahalanobis Distance
                N(i) = (act_meas - meas)' / cov * (act_meas - meas);
                if N(i) < thres
                    % For acceptable error (within the reliability range),
                    % the validity of candidate is set to be true
                    validity(i) = true;
                end
            end
        end
        
        %% Process Validity Array to find acceptable candidate
        function valid_idx = ProcValidity(validity,N)
            num_valid = length(find(validity));
            
            if num_valid == 0 || num_valid > 1
                if num_valid == 0
                    warning('No valid candidate, choosing match from smallest Mahalanobis distance')
                    % Typically, there must be at least one match
                    % There is something wrong if there is no match
                    % * If using INS propagation results directly for
                    % initial states, this may happen.
                    % Using 2-phase optimization is recommended to prevent
                    % situations with 'no match '
                end
                [~,valid_idx] = min(N);
            elseif num_valid == 1
                valid_idx = find(validity);
            end
        end
       
        %% Concatenate normal map segment and extended map segment
        function conc_map = MergeMap(map,ext_map)
            conc_map = struct();
            conc_map.left = [map.left ext_map.left];
            conc_map.right = [map.right ext_map.right];
            conc_map.will_be_matched = map.will_be_matched;
        end

    end
end