function [T] = symt_rotate(W, Q, dsym)
%SYMT_ROTATE Rotates a partially symmetric tensor by a rotation matrix Q
%   [T] = symt_rotate(W, Q)
%
%   Input arguments:
%     W - order-d tensor
%     Q - rotation matrix
%     dsym - number of symmetric dimensions (by default, dsym = d)
%
%   Output arguments:
%     T - the product W *1 Q' *2 Q' ... *dsym Q'; 
  d = ndims(W);
  if (~exist('dsym','var')), dsym = d; end  % default value   
  matr_arg = cell(1, dsym);
  [matr_arg{:}] = deal(Q');
  T = ttm(W, matr_arg, 1:dsym);
end

