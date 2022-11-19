function [lambda_ij] = gradient_lambda_ij_complex(W, i, j, p, d)
%GRADIENT_LAMBDA_ij_COMPLEX Returns a (i,j) gradient element of the cost function
%   [Lambda_ij] = gradient_lambda_ij_complex(W, i, j, p, d)
%
%   Input arguments:
%     W -  complex tensor
%     i,j - pair indeces (need ensure i < j <= n(p))
%     p - mode
%     d - order
%
%   Output argument:
%     lambda_ij - (i,j)-element of the Lambda matrix
%

if (~exist('d','var'))
    d = ndims(W);
end
n     = size(W);
n_min = min(n);
if (p == 1)
    if (j <= n_min)
        lambda_ij = W([i, ones(1, d-1) .* j]) * conj(W([ones(1, d) .* j])) - ...
            conj(W([j, ones(1, d-1) .* i])) * W([ones(1, d) .* i]);
    elseif (i > n_min)
        lambda_ij = 0;
    elseif (i <= n_min) && (j > n_min)
        lambda_ij = - conj(W([j, ones(1, d-1) .* i])) * W([ones(1, d) .* i]);
    end
elseif (p == d)
    if (j <= n_min)
        lambda_ij = W([ones(1, d-1) .* j, i]) * conj(W([ones(1, d) .* j])) - ...
            conj(W([ones(1, d-1) .* i, j])) * W([ones(1, d) .* i]);
    elseif (i > n_min)
        lambda_ij = 0;
    elseif (i <= n_min) && (j > n_min)
        lambda_ij = - conj(W([ones(1, d-1) .* i, j])) * W([ones(1, d) .* i]);
    end
elseif (p > 1) && (p < d)
    if (j <= n_min)
        lambda_ij = W([ones(1, p-1) .* j, i, ones(1, d-p) .* j]) *  conj(W([ones(1, d) .* j])) - ...
            conj(W([ones(1, p-1) .* i, j, ones(1, d-p) .* i])) * W([ones(1, d) .*i]);
    elseif (i > n_min)
        lambda_ij = 0;
    elseif (i <= n_min) && (j > n_min)
        lambda_ij = - conj(W([ones(1, p-1) .* i, j, ones(1, d-p) .* i])) * W([ones(1, d) .* i]);
    end
end
end

