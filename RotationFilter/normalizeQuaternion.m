function x = normalizeQuaternion(x)
% when state is inputted, normalize the quaternion
x(4:7) = x(4:7)/norm(x(4:7));
end