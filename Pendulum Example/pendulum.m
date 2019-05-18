%% clear
clear all;

%%
dt = 0.025;
tend = 60;

tspan = 0:dt:tend;
x0 = [pi/6 0]';

[t,x] = ode45(@(t, x) pendulumState(x), tspan, x0);

%% measurement 
encoder_sigma = .05;
y_meas = x(:, 1) + normrnd(0, encoder_sigma, [numel(t) 1]);

figure
plot(t, y_meas, t, x(:,1))

%% filter
% parameters

L = numel(x(1,:)); % number of states
N = numel(t);

q=1e-2;    %std of process 
r=.1;    %std of measurement
Q=q^2*eye(L); % covariance of process
R=r^2;        % covariance of measurement  

x_hat_k = [pi/12 0.1]'; 
P_k = eye(L);

x_hat = zeros(N, L); % estimate of state

for k=1:(N-1)
    [x_hat_k, P_k] = ukf(@pendulumState, x_hat_k, P_k, @pendulumMeas, y_meas(k+1), Q, R, dt);
    x_hat(k,:) = x_hat_k;
end

figure
plot(t, x_hat(:,1), t, x(:,1))
legend('UKF Estimate', 'Simulated Dynamics')


err_pos = (x_hat(:,1) - x(:,1)).^2;
err_vel = (x_hat(:,2) - x(:,2)).^2;

% figure 
% plot(t, err_pos, t, err_vel)
% legend('Position', 'Velocity')

