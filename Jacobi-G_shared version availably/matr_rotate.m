function [T] = matr_rotate(W, U)
%SYMT_ROTATE Rotates a set of matrices by a rotation matrix U
%   [T] = symt_rotate(W, U)
%
%   Input arguments:
%     W - order-d tensor
%     Q - rotation matrix
%     dsym - number of symmetric dimensions (by default, dsym = d)
%
%   Output arguments:
%     T - the product W *1 U' *2 U.'; 
  n = size(W,1);
  L = size(W,3);
 
  T = zeros(n,n,L);

  for l=1:L
    T(:,:,l) =  U' * W(:,:,l) * U;
  end  
end

