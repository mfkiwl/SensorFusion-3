function R = RightJac(vec)
    mag = norm(vec);
    S = skew(vec);
    if mag < 1e-6
        R = eye(3) - 1/2 * S;
    else
        one_minus_cos = 2 * sin(mag/2) * sin(mag/2);
        R = eye(3) - one_minus_cos/mag^2 * S + (1/mag^2 - sin(mag)/mag^3) * S^2;
    end
end