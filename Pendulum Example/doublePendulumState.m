function dxdt = doublePendulumState(x)
    % pendulum dynamics
    % x = [theta_1 theta_2 theta_1_dot theta_2_dot]
    % http://sophia.dtp.fmph.uniba.sk/~kovacik/doublePendulum.pdf
    
    g = 9.81;
    l1 = 0.1;
    l2 = 0.1;
    m1 = 1;
    m2 = 1;
  
    dxdt = zeros(4,1);
    dxdt(1) = x(3);
    dxdt(2) = x(4);
%     dxdt(3) = -(g/l1).*sin(x(1));  
%     dxdt(4) = -(g/l2).*sin(x(2));  
    dxdt(3) = (-m2*l1*x(4)^2*sin(x(1)-x(2))*cos(x(1)-x(2))+...
                g*m2*sin(x(2))*cos(x(1)-x(2)) - m2*l2*x(4)^2*sin(x(1)...
                - x(2)) - (m1 + m2)*g*sin(x(1)))/(l1*(m1+m2) - ...
                m2*l1*cos(x(1) - x(2))^2);
    dxdt(4) = (-m2*l2*x(4)^2*sin(x(1)-x(2))*cos(x(1)-x(2))+...
                g*sin(x(1))*cos(x(1)-x(2))*(m1 + m2) - ...
                l1*x(4)^2*sin(x(1)-x(2))*(m1 + m2)...
                - (m1 + m2)*g*sin(x(2)))/(l1*(m1+m2) - ...
                m2*l1*cos(x(1) - x(2))^2); 


end
