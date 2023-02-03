function ER = InvMahalanobis(M,cov)
% INVMAHALANOBIS Square root weighting using covariance matrix
% This function is mainly used for weighting in Non-linear Least Squares 
    n = size(cov,1);
    SIG = eye(n)/chol(cov);
    SIG = SIG';
    ER = SIG * M;
end