function y = meas(x)
    % a_p = proper acceleration measured by accelerometer
    % w_m = angular velocity measured by gyroscope
    y = zeros(6,1);
    
    y(1:3) = 9.81;
    y(4:6) = x(11:13);
end