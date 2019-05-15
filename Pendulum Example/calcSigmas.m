function X = calcSigmas(x, P, alpha)
    % x - state
    % P - covariance matrix
    % alpha - initial weighting coefficient 
    
    L = numel(x);
    
    X = zeros(L, 2*L);
    
    lam = eig(P);
    try
        A = alpha*chol(P)';
    catch
        P = P + 1e-5 * eye(L);
        A = alpha*chol(P)';
    end
    Y = repmat(x, 1, L);

    X(:,1:L) = Y+A;
    X(:,(L+1):(2*L)) = Y-A;
    
    
    
%     
%     L = numel(x);
%     
%     X = zeros(L, 2*L + 1);
%     c = sqrt(L/(1-alpha));
%     try
%         A = c*chol(P)';
%     catch
%         P = P + 1e-4 * eye(L);
%         A = c*chol(P)';
%     end
%     Y = repmat(x, 1, L);
% 
%     X(:,1) = x;
%     X(:,2:(L+1)) = Y+A;
%     X(:,(L+2):(2*L + 1)) = Y-A;
end