function R = Exp_map(vec)
    mag = norm(vec);
    S = skew(vec);
    if mag < 1e-6
        R = eye(3) + S;
    else
        one_minus_cos = 2 * sin(mag/2) * sin(mag/2);
        R = eye(3) + sin(mag)/mag * S + one_minus_cos/mag^2 * S^2;
    end
end