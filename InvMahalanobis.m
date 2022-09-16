function ER = InvMahalanobis(M,cov)
    n = size(cov,1);
    SIG = eye(n)/chol(cov);
    SIG = SIG';
    ER = SIG * M;
end