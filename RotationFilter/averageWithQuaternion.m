function xbar = averageWithQuaternion(x, w)
% based off method from:
% Markley, F. Landis, et al. "Averaging quaternions." 
% Journal of Guidance, Control, and Dynamics 30.4 (2007): 1193-1197.
    n = numel(x(:,1));
    
    xbar = zeros(19,1);
    
    xbar([1:3 8:19]) = x([1:3 8:19], :)*w;
    
    Q = x(4:7, :)*diag(sqrt(w));
    A = Q*Q';
    
    [Qbar, Lam] = eig(A);
        
    [~,index] = max(diag(Lam)); 
    qbar = Qbar(:,index); 
    
    xbar(4:7) = qbar;
end