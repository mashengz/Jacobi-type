function [Lambda] = gradient_lambda_complex(W, d)
%GRADIENT_LAMBDA_COMPLEX Returns a Riemannian gradient of the cost function (sum of squares)
%   [Lambda] = gradient_lambda_complex(W)
%
%   Input arguments:
%     W - complex tensor
%     d - order
%
%   Output argument:
%     Lambda - the Lambda matrix, output arguments of type 'cell'
%

if (~exist('d','var'))
    d = ndims(W);
end
n      = size(W);
Lambda = cell(1, d);
for p = 1 : d
    Lambda{p} = zeros(n(p), n(p));
    for i = 1 : (n(p) - 1)
        for j = (i + 1) : n(p)
            Lambda{p}(i, j) = gradient_lambda_ij_complex(W, i, j, p, d);
            Lambda{p}(j, i) = - conj(Lambda{p}(i, j));
        end
    end
end
end

