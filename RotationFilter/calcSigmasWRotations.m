function X = calcSigmasWRotations(x, P, alpha)
    % calculates the sigma points with multiplicative quaternion update for
    %   q = x(4:7)
    % x - mean estimate of the state
    % P - covariance matrix
    % alpha - weighting coefficient (higher alpha spreads out distribution)
    
    L = numel(x);
    
    X = zeros(L, 2*(L-1));
    
%     lam = eig(P);
    
%     A = alpha*chol(P)';
    try
        A = alpha*chol(P)';
%         disp(cond(P));
    catch
        % numerical conditioning 
        P = P + 1e-5 * eye(L-1);
        A = alpha*chol(P)';
    end

    for i=1:(L-1)       
        % do simple additive sigma point calculation
        X([1:3 8:19], i) = x([1:3 8:19]) + A([1:3 7:18], i);
        X([1:3 8:19], i+(L-1)) = x([1:3 8:19]) - A([1:3 7:18], i);
        
        % do quaternion multiplicative update
        phi = A(4:6,i);
        dq_plus = quat_exp(phi);
        dq_neg = quat_exp(-phi);
        
        X(4:7, i) = quat_prod(x(4:7), dq_plus);
        X(4:7, i+(L-1)) = quat_prod(x(4:7), dq_neg);
    end
   
end