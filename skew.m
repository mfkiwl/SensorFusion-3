function S = skew(vec)
    % Returns skew symmetric matrix
    if length(vec) ~= 3
        error('3*1 vector form is expected')
    end
    S = [0 -vec(3) vec(2);
         vec(3) 0 -vec(1);
         -vec(2) vec(1) 0];
end