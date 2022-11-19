function [Lambda_p] = gradient_hess_h_complex_r(A, M, Q, p, d)
%GRADIENT_Hess_H_COMPLEX Returns a (i,j) gradient and hess element of the cost function
%   [Lambda_p] = gradient_hess_h_complex_r(A, M, Q, p, d)
%
%   Input arguments:
%     A -  complex tensor
%     M,Q - matrices (in fact, n(p) \times n(p))
%     p - mode
%     d - order
%
%   Output arguments:
%     Lambda_p - the mode-p Rie-gradient Lambda
%

if (~exist('d','var'))
    d = ndims(A);
end
n = size(A);
n_min = min(n);
MQ = cellfun(@mtimes, Q, M, 'un', 0);
W = tenmat_prod(A, MQ, d);
v = zeros(n_min, 1);
matr_arg      = cell(1, d);
[matr_arg{:}] = deal(':');
for i = 1 : n_min % check it, n_min -> r(p) ???
    [matr_arg{1:d}] = deal(i);
    v(i) = W(matr_arg{:});
end
D = tendiag(v,n);
W_p = tens2mat(D.data, p);
B = tenmat_prod_p(A, MQ, d, p);
B_p = tens2mat(B.data, p);
% compute Lambda matrix
Lambda_p = Q{p}'*(B_p*W_p')*M{p};
Lambda_p = Lambda_p - Lambda_p';
end