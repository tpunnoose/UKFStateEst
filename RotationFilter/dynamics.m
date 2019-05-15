function xdot = dynamics(x)
    % rotating cube
    m = 1.0;
    s = 1;

    xdot = zeros(13,1);
    
    
    J1 = 1;
    J2 = 2;
    J3 = 3;
    
    J = diag([J1 J2 J3]);
    
    q = x(4:7);
    v = x(8:10);
    om = x(11:13);
    
    % velocity
    xdot(1:3) = v;
    
    % quaternion kinematics
    om_hat = [0 om']';
    xdot(4:7) = 0.5*quat_prod(q,om_hat);
    
    % acceleration
    xdot(8:10) = 0; % no external forces
    
    % angular acceleration
    xdot(11:13) = -inv(J)*cross(om, J*om); % no external torques
end