%% Average Test

clear all;

N = 1000;

q0 = [1 0 0 0]'

x0 = [zeros(3,1); q0; zeros(6,1)];

Q = repmat(x0, 1, N);

r = 10;

V = r^2*randn(13, N);

M = Q+V;

w = 1/N*zeros(N,1);

xbar = averageWithQuaternion(M,w);

qbar = xbar(4:7)

