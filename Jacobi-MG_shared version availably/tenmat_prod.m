function [T] = tenmat_prod(W, Q, d)
% TENMAT_PROD Returns a tensor by some transform matrices Q{1:d}
%    [T] = tenmat_prod(W, Q, d)
%
%   Input arguments:
%     W - tensor
%     Q - some transform matrice (the argument must be cell matrix)
%     d - order
%
%   Output argument:
%     T - the product W *1 Q1' *2 Q2' ... *d Qd'
%

if (~exist('d','var'))
    d = ndims(W);
end
mat_arg = cell(1, d);
for i = 1 : d
    [mat_arg{i}] = deal(Q{i}');
end
T = ttm(W, mat_arg, 1:d);
end

