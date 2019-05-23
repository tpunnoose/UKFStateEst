%% clear
clear all;

%% Simulate Dynamics
dt = 0.025;
tend = 60;

tspan = 0:dt:tend;

r0 = [0 0 0]';
q0 = [1 0 0 0]';
v0 = [0 0.1 0]';
om0 = [1 1 1]';
bf0 = [0 0 0]'; % accelerometer bias
bw0 = [0 0 0]'; % gyro bias

x0 = [r0; q0; v0; om0; bf0; bw0];

N = tend/dt + 1;
X = zeros(19, N-1);
X(:,1) = x0;

for i=2:(N-1)
    k1 = dt*dynamics_sim(X(:,i-1));
    k2 = dt*dynamics_sim(X(:,i-1)+k1/2);
    k3 = dt*dynamics_sim(X(:,i-1)+k2/2);
    k4 = dt*dynamics_sim(X(:,i-1)+k3);
    X(:,i) = X(:,i-1) + (k1+2*k2+2*k3+k4)/6;
    
    % renormalize quaternion
    X(4:7,i) = X(4:7,i)/norm(X(4:7,i));        
end

%% Sensor Data
Y = zeros(6, N);

sigma_accel = 0.05;
sigma_gyro = 0.1;
for i=1:(N-1)
   Y(:, i) = meas(X(:,i));
   Y(1:3, i) = Y(1:3, i) + normrnd(0, sigma_accel, 3, 1);
   Y(4:6, i) = Y(4:6, i) + normrnd(0, sigma_gyro, 3, 1);
end

%% Filtering
L = numel(X(:,1)); % number of states
M = numel(Y(:,1)); % number of measurements
N = numel(tspan);

q=0.01;    %std of process 
r=0.01;    %std of measurement
Q=q^2*eye(L-1); % covariance of process
R=r^2*eye(M);        % covariance of measurement  

x_hat_k = x0;
x_hat_k(4:7) = quat_exp(pi/12 * [1 0 0]');

% because the quaternion has 3 DOF and 4 elements
P_k = .1*eye(L-1);

x_hat = zeros(N-1, L); % estimate of state

for k=1:(N-1)
    [x_hat_k, P_k] = mukf(@dynamics_est, x_hat_k, P_k, ...
                        @meas, Y(:, k+1), Q, R, dt);
    x_hat(k,:) = x_hat_k;
end


%% Find Quaternion Error

mag = zeros(N-1,1);

for i=1:(N-1)
    dq = quat_prod((X(4:7, i).*[1 -1 -1 -1]'), x_hat(i, 4:7)');
    
    mag(i) = abs(wrapToPi(norm(quat_log(dq))));
end

figure
plot(tspan(1:(end-1)), mag)
%% Plot

figure
plot(tspan(1:(end-1)), X(4,:), tspan(1:(end-1)), x_hat(:,4))
