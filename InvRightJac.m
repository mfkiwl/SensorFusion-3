function R = InvRightJac(vec)
    mag = norm(vec);
    S = skew(vec);
    if mag < 1e-6
        R = eye(3);
    else
        R = eye(3) + 1/2 * S + (1/mag^2 + (1+cos(mag))/(2*mag*sin(mag))) * S^2;
    end
end