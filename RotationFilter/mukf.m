function [x, P] = mukf(f, x_prev, P_prev, h, z_meas, Q, R, dt)
    % f - state function (CT) -> integrated by RK4
    % x_prev - previous state estimate 
    % P_prev - previous covariance
    % h - measurement function (DT)
    % z_meas - new measurement
    % Q - process covariance
    % M - measurement covariance
    % dt - time step
    L = numel(x_prev); % number of elements in state vector

    alpha = .1;
        
    % generate sigma points:
    X_x = calcSigmasWRotations(x_prev, P_prev, alpha);
    
    n = size(X_x,2); % number of sigma points
    x_p = zeros(L,n);
    
    
    % weighting vectors:
    w_m = 1/(n)*ones(1, n)';
    w_c = 1/alpha^2 * w_m;
    
    % propogate sigma points through state function and integrate via
    % midpoint integration:
%     for k=1:n 
%         dxdt_1 = f(X_x(:,k)); 
%         dxdt_2 = f(X_x(:,k) + dt/2*dxdt_1);
%         x_p(:,k) = X_x(:,k) + dt*dxdt_2;
%     end

    %rk4 integration:
    for k=1:n 
        k1 = dt*f(X_x(:,k));
        k2 = dt*f(X_x(:,k)+k1/2);
        k3 = dt*f(X_x(:,k)+k2/2);
        k4 = dt*f(X_x(:,k)+k3);
        x_p(:,k) = X_x(:,k) + (k1+2*k2+2*k3+k4)/6;
        x_p(:,k) = normalizeQuaternion(x_p(:,k));
    end

    xbar = averageWithQuaternion(x_p, w_m);
    
    z_p = zeros(numel(z_meas), n);
    % run propogated points through measurement model 
    for k=1:n
        z_p(:,k) = h(x_p(:,k));
    end
    
    zbar = zeros(numel(z_meas), 1);
    % calculate predicted measurement
    for k=1:n
        zbar = zbar + w_m(k)*z_p(:,k);
    end
    
    % calculate covariances
    
    dx_p = calcStateDiff(x_p, xbar);
    
    % calculate covariance of predicted state
    Pbar = Q;
    for j = 1:n
        Pbar = Pbar + w_c(j)*(dx_p(:,j))*(dx_p(:,j))';
    end
    
    % calculate covariance of predicted measurement
    Pzz = zeros(numel(z_meas),numel(z_meas));
    for j = 1:n
        Pzz = Pzz + w_c(j)*(z_p(:,j) - zbar)*(z_p(:,j) - zbar)';
    end
    
    % calculate cross-covariance
    Pxz = zeros(L-1,numel(z_meas));
    for j = 1:n
        Pxz = Pxz + w_c(j)*(dx_p(:,j))*(z_p(:,j) - zbar)';
    end
    
    %Innovation
    nu = z_meas - zbar; 
    
    %Innovation Covariance
    S = Pzz + R; 
    
    %Kalman Gain
    K = Pxz/S; 
    
    x_up = K*nu;
    
    % multiplicative quaternion update
    x([1:3 8:19]) = (xbar([1:3 8:19]) + x_up([1:3 7:18]));
    x(4:7) = quat_prod(xbar(4:7), quat_exp(x_up(4:6)));
    
    x = x';
    
    P = Pbar - K*S*K';
end