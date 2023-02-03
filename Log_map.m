function vec = Log_map(R)
    tr_3 = trace(R) - 3;

    if tr_3 < -1e-6
        psi = acos(1/2 * (trace(R) - 1));
        magnitude = psi / (2*sin(psi));
    else
        magnitude = 0.5 - 1/12 * tr_3 + tr_3^2 /60;
    end

    vec = magnitude * [R(3,2) - R(2,3);R(1,3) - R(3,1);R(2,1) - R(1,2)];
end