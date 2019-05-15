function [xhat, Phat, t, Khist] = UKF_zac(x0, P0, Q, R, zhist, dt)
    Nx = length(x0);

    %Tuning parameters for sigma points
    alpha = .05;
    beta = 2; %2 is optimal for a Gaussian distribution
    kappa = -3; %choose kappa so that kappa + Nx = 3

    lambda = alpha^2*(kappa + Nx) - Nx;

    %Compute sigma point weights
    w_m = 1/(2*(Nx+lambda))*ones(2*Nx+1,1);
    w_m(1) = lambda/(Nx+lambda);
    w_c = w_m;
    w_c(1) = lambda/(Nx+lambda) + 1 - alpha^2 + beta^2;

    %Initialize arrays
    t = (0:(length(zhist)-2))*dt;
    xhat = zeros(2,length(t)); 
    Phat = zeros(2,2,length(t));
    Khist = zeros(2,1,length(t));
    
    xhat(:,1) = x0;
    Phat(:,:,1) = P0;

    for k = 1:length(t)

        %Generate sigma points
        S = chol(Phat(:,:,k))';
        chi = zeros(Nx,2*Nx+1);
        chi(:,1) = x0;
        for j = 1:Nx
            chi(:,1+j) = chi(:,1) + sqrt(Nx+lambda)*S(:,j);
            chi(:,1+Nx+j) = chi(:,1) - sqrt(Nx+lambda)*S(:,j);
        end

        %Propagate sigma points forward one time step
        x_p = zeros(Nx,1+2*Nx);
        for j = 1:(1+2*Nx)
            x_p(:,j) = pendulumState(chi(1:Nx,j));
        end

        %Run propagated points through measurement model
        z_p = zeros(size(zhist,1), 2*Nx+1);
        for j = 1:(2*Nx+1)
            z_p(:,j) = pendulumMeas(x_p(:,j));
        end

        %Calculate predicted state
        xbar = zeros(Nx,1);
        for j = 1:length(w_m)
            xbar = xbar + w_m(j)*x_p(:,j);
        end

        %Calculate predicted measurement
        zbar = zeros(size(zhist,1),1);
        for j = 1:length(w_m)
            zbar = zbar + w_m(j)*z_p(:,j);
        end

        %Calculate predicted covariances
        Pbar = Q;
        for j = 1:length(w_c)
            Pbar = Pbar + w_c(j)*(x_p(:,j) - xbar)*(x_p(:,j) - xbar)';
        end
        Pzz = zeros(size(zhist,1),size(zhist,1));
        for j = 1:length(w_c)
            Pzz = Pzz + w_c(j)*(z_p(:,j) - zbar)*(z_p(:,j) - zbar)';
        end
        Pxz = zeros(Nx,size(zhist,1));
        for j = 1:length(w_c)
            Pxz = Pxz + w_c(j)*(x_p(:,j) - xbar)*(z_p(:,j) - zbar)';
        end

        %Measurement Update
        nu = zhist(:,k+1) - zbar; %Innovation
        S = Pzz + R; %Innovation Covariance
        K = Pxz/S; %Kalman Gain
        Khist(:,:,k) = K;

        xhat(:,k+1) = xbar + K*nu;
        Phat(:,:,k+1) = Pbar - K*S*K';    

    end
end

