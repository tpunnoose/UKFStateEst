function dxdt = pendulumState(x)
    % pendulum dynamics
    g = 9.81;
    l = 0.1;
    
    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = -(g/l).*sin(x(1));    
end
