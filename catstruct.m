function newS = catstruct(S)
% CATSTRUCT(S) assumes that all the structures inside the cell S have
% identical fields and concatenate them.
% This script is not for generalized usage. 
%
% Implemented by JinHwan Jeon, 2023

    n = length(S); % Number of structures given
    ref = S{1}; % Pick reference structure to sample fieldnames
    fields = fieldnames(ref);
    
    newS = struct();
    for i=1:numel(fields)
        field = fields{i}; % rtvel, rtgps, ...
        fieldX2s = fieldnames(ref.(field));        
        newS.(field) = struct();
        M = zeros(numel(fieldX2s),1);
        N = zeros(numel(fieldX2s),1);
        
        for j=1:numel(fieldX2s)
            fieldX2 = fieldX2s{j}; % t, lat, lon, ...
            M(j) = size(ref.(field).(fieldX2),1);
            N(j) = 0;
            for k=1:n
                sampleS = S{k};
                N(j) = N(j) + size(sampleS.(field).(fieldX2),2);
            end
        end
        
        for j=1:numel(fieldX2s)
            fieldX2 = fieldX2s{j}; % t, lat, lon, ...
            newS.(field).(fieldX2) = zeros(M(j),N(j));
            cnt = 0;
            for k=1:n
                sampleS = S{k};
                N_ = size(sampleS.(field).(fieldX2),2);
                newS.(field).(fieldX2)(:,cnt+1:cnt+N_) = sampleS.(field).(fieldX2);
                cnt = cnt + N_;
            end
        end        
    end
end