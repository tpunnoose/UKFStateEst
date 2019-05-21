function ql = quat_log(q)
% returns the quaternion log of the input quaternion
% via 2.4.4 of https://www.iri.upc.edu/people/jsola/JoanSola/objectes/notes/kinematics.pdf
eps = 0.01;

if(norm(q(2:4)) > eps)
    phi = 2*acos(q(1));%2*atan2(norm(q(2:4)), q(1));
    u = q(2:4)/norm(q(2:4));
    
    ql = phi*u;
else % use taylor series approximation so small angles don't blow up
    ql = 2*q(2:4)/q(1)*(1 - norm(q(2:4))^2/(3*q(1)^2));
end
end