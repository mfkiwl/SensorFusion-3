function [I, J, V] = sparseFormat(rows, cols, values)
    m = length(rows); n = length(cols);
    if ne(size(values), [m n])
        error('Value matrix format does not match row, column information')
    end
    
    I = zeros(1,m*n); J = zeros(1,m*n); V = zeros(1,m*n);
    for i=1:m
        for j=1:n
            I((i-1)*n+j) = rows(i);
            J((i-1)*n+j) = cols(j);
            V((i-1)*n+j) = values(i,j);
        end
    end
end