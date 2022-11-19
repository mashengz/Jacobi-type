function [Lambda] = gradient_lambda_complex_r(A, M, Q, d)
%GRADIENT_LAMBDA_COMPLEX Returns a Riemannian gradient of the cost function (sum of squares)
%   [Lambda] = gradient_lambda_complex_r(W)
%
%   Input arguments:
%     A - complex tensor
%     M,Q - matrices (in fact, n(p) \times n(p))
%     d - order
%
%   Output argument:
%     Lambda - the Lambda matrix, output arguments of type 'cell'
%

if (~exist('d','var'))
    d = ndims(A);
end
Lambda = cell(1, d);
for p = 1 : d
    Lambda{p} = gradient_hess_h_complex_r(A, M, Q, p, d);
end
end

