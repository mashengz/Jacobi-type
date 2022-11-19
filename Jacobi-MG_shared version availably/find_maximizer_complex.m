function [x, f] = find_maximizer_complex(W, i, j, p, d)
%FIND_MAXIMIZER_COMPLEX Returns an optimal solution
%   [x, f] = find_maximizer_complex(W, i, j, p, d)
%
%   Input arguments:
%     W -  complex tensor
%     i,j - pair indeces
%     p - mode
%     d - order
%
%   Output arguments:
%     x - an optimal solution satisfies x(1) >= 0
%     f - increment, h(x) - h(0)
%

if (~exist('d','var'))
    d = ndims(W);
end
[lambda_ij] = gradient_lambda_ij_complex(W, i, j, p, d);
n           = size(W);
n_min       = min(n);
if (p == 1)
    if (j <= n_min)
        phi_1 = abs(W([ones(1, d) .* i]))^2 + abs(W([ones(1, d) .* j]))^2;
        phi_2 = abs(W([j, ones(1, d-1) .* i]))^2 +...
            abs(W([i, ones(1, d-1) .* j]))^2;
    elseif (i > n_min)
        phi_1 = 0;
        phi_2 = 0;
    elseif (i <= n_min) && (j > n_min)
        phi_1 = abs(W([ones(1, d) .* i]))^2;
        phi_2 = abs(W([j, ones(1, d-1) .* i]))^2;
    end
elseif (p == d)
    if (j <= n_min)
        phi_1 = abs(W([ones(1, d) .* i]))^2 + abs(W([ones(1, d) .* j]))^2;
        phi_2 = abs(W([ones(1, d-1) .* i, j]))^2 +...
            abs(W([ones(1, d-1) .* j, i]))^2;
    elseif (i > n_min)
        phi_1 = 0;
        phi_2 = 0;
    elseif (i <= n_min) && (j > n_min)
        phi_1 = abs(W([ones(1, d) .* i]))^2;
        phi_2 = abs(W([ones(1, d-1) .* i, j]))^2;
    end
elseif (p > 1) && (p < d)
    if (j <= n_min)
        phi_1 = abs(W([ones(1, d) .* i]))^2 + abs(W([ones(1, d) .* j]))^2;
        phi_2 = abs(W([ones(1, p-1) .* i, j, ones(1, d-p) .* i]))^2 +...
            abs(W([ones(1, p-1) .* j, i, ones(1, d-p) .* j]))^2;
    elseif (i > n_min)
        phi_1 = 0;
        phi_2 = 0;
    elseif (i <= n_min) && (j > n_min)
        phi_1 = abs(W([ones(1, d) .* i]))^2;
        phi_2 = abs(W([ones(1, p-1) .* i, j, ones(1, d-p) .* i]))^2;
    end
end
phi_3 = - real(lambda_ij);
phi_4 = - imag(lambda_ij);

% Create M
M = [phi_1 - phi_2, phi_3, phi_4; phi_3, 0, 0; phi_4, 0, 0];

% Find the largest eigenvalue and eigenvector
[u, s] = eig(M);
[maxeigv, ind] = max(sum(s));
% Or use an explicit expression
% t = sqrt(M(1,1)^2 + 4*(phi_3^2 + phi_4^2));
% maxeigv = (M(1,1) + t)/2;

% Calculate increment
f = maxeigv - M(1,1);

% Ensure the output solution satisfies x(1) >= 0
x = u(:, ind);

% Or use an explicit expression
% if (phi_3 == 0) && (phi_4 == 0)
%     x = [1,0,0]';
% elseif (phi_3 ~= 0) && (phi_4 == 0)
%     y = [maxeigv/phi_3,1,0]';
%     x = y/norm(y);
% elseif (phi_3 ~= 0) && (phi_4 ~= 0)
%     y = [maxeigv/phi_4,phi_3/phi_4,1]';
%     x = y/norm(y);
% end
if (x(1) < 0)
    x = - x;
end
end
