%% clear
clear all;

%%
dt = 0.01;
tend = 60;

tspan = 0:dt:tend;
x0 = [pi/12 pi/12 0 0]';

[t,x] = ode45(@(t, x) doublePendulumState(x), tspan, x0);

%% measurement 
encoder_sigma = 0.05;
y_meas = x(:, 1:2) + normrnd(0, encoder_sigma, [numel(t) 2]);

figure
plot(t, y_meas(:,1), t, x(:,1))

%% filter
% parameters

L = numel(x(1,:)); % number of states
N = numel(t);

q=1e-3;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(L); % covariance of process
R=r^2*eye(L/2);        % covariance of measurement  

x_hat_k = [pi/6 pi/6 0 0]'; 
P_k = eye(L);

x_hat = zeros(N, L); % estimate of state

for k=1:(N-1)
    [x_hat_k, P_k] = ukf(@doublePendulumState, x_hat_k, P_k, ...
                        @doublePendulumMeas, y_meas(k+1, :)', Q, R, dt);
    x_hat(k,:) = x_hat_k;
end

figure
plot(t, x_hat(:,1), t, x(:,1))
legend('UKF Estimate', 'Simulated Dynamics')


err = (x_hat(:,1) - x(:,1)).^2;
figure 
plot(t, err)


