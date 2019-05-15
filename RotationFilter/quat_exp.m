function q = quat_exp(phi)
% turn rotation vector into quaternion 
    q = zeros(4,1);
    mag = norm(phi);
    eps = 0.001;

    if mag > eps
        q(1) = cos(0.5*mag);
        q(2:4) = sin(0.5*mag) * phi/mag;
    else
        q(1) = 1;
    end
end