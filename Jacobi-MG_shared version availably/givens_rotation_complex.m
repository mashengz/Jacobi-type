function [G] = givens_rotation_complex(n, i, j, p, x)
%GIVENS_ROTATION_COMPLEX Creates a complex Givens rotation matrix
%   [G] = givens_rotation_complex(n, i, j, p, x)
%
%   Input arguments:
%     n - size of the complex tensor
%     i,j - pair indeces
%     p - mode
%     x - an optimal solution of the rotation
%
%   Output argument:
%     G - a complex Givens matrix
%

G  = eye(n(p));
c  = x(1);
s1 = x(2);
s2 = x(3);
G([i j],[i,j]) = [c, -(s1 + 1i*s2); (s1 - 1i*s2)  c];
end

