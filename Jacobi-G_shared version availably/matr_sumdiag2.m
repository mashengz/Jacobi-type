function [f] = matr_sumdiag2(W, dsym)
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
  L = size(W, 3);
  n = size(W, 1);
 
  f = 0;
  for l=1:L
    f = f + norm(squeeze(diag(W(:,:,l))))^2;
  end  
end

