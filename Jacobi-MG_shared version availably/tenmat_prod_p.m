function [T] = tenmat_prod_p(W, Q, d, p)
% TENMAT_PROD_P Returns a tensor by some transform matrices Q{1:p-1;p+1:d}
%    [T] = tenmat_prod_p(W, Q, d)
%
%   Input arguments:
%     W - tensor
%     Q - some transform matrice (the argument must be cell matrix)
%     d - order
%
%   Output argument:
%     T - the product W *1 Q1'... *(p-1) Q(p-1)'*(p+1) Q(p+1)' ... *d Qd'
%


if (~exist('d','var'))
    d = ndims(W);
end
mat_arg = cell(1, d-1);
if (p == 1)
    for i = 2 : d
        [mat_arg{i}] = deal(Q{i}');
    end
    T = ttm(W, mat_arg, 2:d);
elseif (p == d)
    for i = 1 : d-1
        [mat_arg{i}] = deal(Q{i}');
    end
    T = ttm(W, mat_arg, 1:d-1);
else
    mat_arg_1 = cell(1, p-1);
    mat_arg_2 = cell(1, d-p);
    for i = 1 : p-1
        [mat_arg_1{i}] = deal(Q{i}');
    end
    C = ttm(W, mat_arg_1, 1:p-1);
    for i = p+1 : d
        [mat_arg_2{i}] = deal(Q{i}');
    end
    T = ttm(C, mat_arg_2, p+1:d);
end
end

