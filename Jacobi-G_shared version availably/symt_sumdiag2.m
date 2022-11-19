function [f] = symt_sumdiag2(W, dsym)
%SYMT_ROTATE Calculates the sum of squared elements on the diagonals of the
%symmetric slices.
%   [f] = symt_sumdiag2(W)
%
%   Input arguments:
%     W - order-d  tensor
%     dsym - number of symmetric dimensions (by default, dsym = d)
%
%   Output arguments:
%      f = ||W(1,...,1,:)||^2 + ... + ||W(n,...,n,:)||^2; 
  d = ndims(W);
  n = size(W, 1);
  if (~exist('dsym','var')), dsym = d; end  % default value   
  
  matr_arg = cell(1, d);
  [matr_arg{:}] = deal(':');
 
  f = 0;
  for i=1:n
    [matr_arg{1:dsym}] = deal(i);
    f = f + norm(W(matr_arg{:}))^2;
  end  
end

