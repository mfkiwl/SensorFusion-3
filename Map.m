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
    
    properties (Access = public)
        states
        lane
        segments = {}
        lane_prob_thres
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
                for j = lb:ub-1
                    R = obj.states{j}.R;
                    p = obj.states{j}.p;
                    map.left = [map.left p + R * obj.states{j}.left(:,1)];
                    map.right = [map.right p + R * obj.states{j}.right(:,1)];
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
                
                % For last state of the segment, add all previewed points
                R = obj.states{ub}.R;
                p = obj.states{ub}.p;
                map.left = [map.left p + R * obj.states{ub}.left];
                map.right = [map.right p + R * obj.states{ub}.right];

                obj.segments = [obj.segments {map}];
            end
        end
        
        %% Lane Point data association (within segment)
        % * Any way to make this more efficient?
        % 
        function obj = associate(obj)
            n = size(obj.lane.FactorValidIntvs,1);
            for i=1:n
                map = obj.segments{i};
                lb = obj.lane.FactorValidIntvs(i,1);
                ub = obj.lane.FactorValidIntvs(i,2); 
                assocL = zeros(ub - lb,obj.lane.prev_num-1);
                assocR = assocL;
                

                for j=lb:ub-1
                    if map.will_be_matched(j)
                        R = obj.states{j}.R;
                        P = obj.states{j}.P;
                        LEFT = obj.states{j}.left;
                        RIGHT = obj.states{j}.right;
                        % Perform frame transformation w.r.t vehicle frame
                        map_b = frameTransformation(R,P,map);  
                        xL = map_b.left(1,:);
                        xR = map_b.right(1,:);
                        for k=2:obj.lane.prev_num
                            idxL = find(xL > 10*(k-1) - 5 & xL < 10*(k-1) + 5); % may need to change value 5 to other valid values
                            idxR = find(xR > 10*(k-1) - 5 & xR < 10*(k-1) + 5);
                            
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

                                    Rt = obj.states{obj.getStateIdx(idxP,i)}.R;
                                    Pt = obj.states{obj.getStateIdx(idxP,i)}.P;

                                    lane_idxP = obj.lane.state_idxs(obj.getStateIdx(idxP,i));
                                    lane_idxN = obj.lane.staet_idxs(obj.getStateIdx(idxN,i));
                                    x_ = [xL(idxP) xL(idxN)];
                                    stdy_ = [obj.lane.lystd(lane_idxP,k) obj.lane.lystd(lane_idxN,k)];
                                    stdz_ = [obj.lane.lzstd(lane_idxP,k) obj.lane.lystd(lane_idxN,k)];
                                    y_ = [obj.lane.ly(lane_idxP,k) obj.lane.ly(lane_idxN,k)];
                                    z_ = [obj.lane.lz(lane_idxP,k) obj.lane.lz(lane_idxN,k)];

                                    stdy_inter = interp1(x_,stdy_,10*(k-1));
                                    stdz_inter = interp1(x_,stdz_,10*(k-1));
                                    y_inter = interp1(x_,y_,10*(k-1));
                                    z_inter = interp1(x_,z_,10*(k-1));
                                    meas_inter = [y_inter;z_inter];
                                    cov_inter = diag([stdy_inter^2,stdz_inter^2]);

                                    cand = struct();
                                    cand.meas = meas_inter;
                                    cand.cov = cov_inter;
                                    % Need to convert actual world frame
                                    % lane point coords to matching frame

                                    cand.actual_meas = Rt' * (P + R * LEFT(:,k) - Pt);
                                    candL = [candL {cand}];
                                end                                
                            end

                            % Right

                            for t=1:length(idxR)-1
                                idxP = idxR(t); idxN = idxR(t+1);
                                if xR(idxP) < 10*(k-1) && xR(idxN) >= 10*(k-1)
                                    % If match candidate is found,
                                    % interpolate and find resample point.
                                    % There may be multiple candidates
                                    
                                    Rt = obj.states{obj.getStateIdx(idxP,i)}.R;
                                    Pt = obj.states{obj.getStateIdx(idxP,i)}.P;

                                    lane_idxP = obj.lane.state_idxs(obj.getStateIdx(idxP,i));
                                    lane_idxN = obj.lane.staet_idxs(obj.getStateIdx(idxN,i));
                                    x_ = [xR(idxP) xR(idxN)];
                                    stdy_ = [obj.lane.rystd(lane_idxP,k) obj.lane.rystd(lane_idxN,k)];
                                    stdz_ = [obj.lane.rzstd(lane_idxP,k) obj.lane.rystd(lane_idxN,k)];
                                    y_ = [obj.lane.ry(lane_idxP,k) obj.lane.ry(lane_idxN,k)];
                                    z_ = [obj.lane.rz(lane_idxP,k) obj.lane.rz(lane_idxN,k)];

                                    stdy_inter = interp1(x_,stdy_,10*(k-1));
                                    stdz_inter = interp1(x_,stdz_,10*(k-1));
                                    y_inter = interp1(x_,y_,10*(k-1));
                                    z_inter = interp1(x_,z_,10*(k-1));
                                    meas_inter = [y_inter;z_inter];
                                    cov_inter = diag([stdy_inter^2,stdz_inter^2]);

                                    cand = struct();
                                    cand.meas = meas_inter;
                                    cand.cov = cov_inter;
                                    cand.actual_meas = Rt' * (P + R * RIGHT(:,k) - Pt);
                                    candR = [candR {cand}];
                                end                                
                            end

                            % Verify candidates by performing Chi-Square
                            % Test. 
                            
                        end
                    end
                end

                % For last index
            end

        end
        
        %% Get State number from map index and segment number
        % Except for last state indices
        function state_idx = getStateIdx(obj,map_idx,seg_idx)
            curr_seg_intv = obj.lane.FactorValidIntvs(seg_idx,1:2);
            bds = curr_seg_intv(1):curr_seg_intv(2);
            state_idx = bds(map_idx);
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
        function validity = ChiSquareTest(cand)
            n = length(cand);
            validity = false(1,n);
            
            % 99% Reliability for 2 DOF system
            thres = chi2inv(0.99,2); 

            for i=1:n
                act_meas = cand{i}.actual_meas(2:3);
                meas = cand{i}.meas;
                cov = cand{i}.cov;

                N = (act_meas - meas) / cov * (act_meas - meas);
                if N < thres
                    % For acceptable error (within the reliability range),
                    % the validity of candidate is set to be true
                    validity(i) = true;
                end
            end
        end

    end
end