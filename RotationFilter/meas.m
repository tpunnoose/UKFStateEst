function y = meas(x)
    % a_p = proper acceleration measured by accelerometer
    % w_m = angular velocity measured by gyroscope
    y = zeros(6,1);
    
    q_g_w = [0 0 0 9.81]';
    
    q_g_b = quat_prod(x(4:7), quat_prod(q_g_w, x(4:7).*[1 -1 -1 -1]));
    
    y(1:3) = q_g_b(2:4);
    y(4:6) = x(11:13);
end