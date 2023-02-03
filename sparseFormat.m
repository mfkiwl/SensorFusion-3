function [I, J, V] = sparseFormat(rows, cols, values)
% sparseFormat returns the matrix indices and corresponding value array in  
% sparse matrix format. This function is mainly used for constructing
% sparse Jacobian matrix efficiently.
%
% [I, J, V] = sparseFormat(rows,cols,values)
% A = sparse(I,J,V,m,n) --> sparse matrix can be constructed directly

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