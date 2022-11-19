function [gradient_h0] = gradient_hess_h_complex(W, i, j, p, d)
%GRADIENT_Hess_H_COMPLEX Returns a (i,j) gradient and hess element of the cost function
%   [gradient_h0] = gradient_hess_h_complex(W, i, j, d, p)
%
%   Input arguments:
%     W -  complex tensor
%     i,j - pair indeces (in fact, i < j <= n(p))
%     p - mode
%     d - order
%
%   Output arguments:
%     gradient_h0 - a (i,j)-element of the mode-p h_gradient
%

if (~exist('d','var'))
    d = ndims(W);
end
Lambda = gradient_lambda_complex(W, d);
gradient_h0 = sqrt(2)*abs(Lambda{p}(i,j));
end