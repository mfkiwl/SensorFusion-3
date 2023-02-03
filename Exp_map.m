function R = Exp_map(vec)
% Exp_map Exponential Mapping for R3 vector to SO(3) Rotational Matrix
% R = Exp_map(vec)
%
% This function is used to compute incremental rotation using 3D vectors in
% IMU Preintegration and various optimization residuals/jacobians.

    mag = norm(vec);
    S = skew(vec);
    if mag < 1e-6
        R = eye(3) + S;
    else
        one_minus_cos = 2 * sin(mag/2) * sin(mag/2);
        R = eye(3) + sin(mag)/mag * S + one_minus_cos/mag^2 * S^2;
    end
end