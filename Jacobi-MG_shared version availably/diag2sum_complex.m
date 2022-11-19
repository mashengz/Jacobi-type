function [f] = diag2sum_complex(W, d)
%DIAG2SUM_COMPLEX Calculates the sum of squared elements on the diagonals of the
% smallest slices
%   [f] = diag2sum_complex(W, d)
%
%   Input arguments:
%     W - complex tensor
%     d - order
%
%   Output argument:
%      f - |W(1,...,1)|^2 + ... + |W(n,...,n)|^2
%

if (~exist('d','var'))
    d = ndims(W);
end
n             = min(size(W));
matr_arg      = cell(1, d);
[matr_arg{:}] = deal(':');
f             = 0;
for i = 1 : n
    [matr_arg{1:d}] = deal(i);
    f               = f + abs(W(matr_arg{:}))^2;
end
end

