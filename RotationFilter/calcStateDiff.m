function dx = calcStateDiff(x, xbar)
    L = numel(xbar);
    N = numel(x(1,:));
    
    dx = zeros(L-1, N);
    q0_inv = xbar(4:7).*[1 -1 -1 -1]';
    for i=1:N
        dx([1:3 7:18], i) = x([1:3 8:19], i) - xbar([1:3 8:19]);
        dq_i = quat_prod(q0_inv, x(4:7, i));
        dx(4:6, i) = quat_log(dq_i);
    end
end