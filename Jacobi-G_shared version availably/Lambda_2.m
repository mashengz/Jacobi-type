function [Lambda] = Lambda_2(W)
%LAMBDA_2 Returns the matrix Lambda for joint diafonalization
%   [Lambda] = Lambda_2(W)
%
%   Input arguments:
%     W - nxnxL tensor
%
%   Output arguments:
%     Lambda - the Lambda matrix
  n = size(W,1);
  L = size(W,3);
  
  Lambda = zeros(n,n);
  
  for l=1:L
    Wl = W(:,:,l);
    Wdiag = conj(diag(Wl));
    Lambda = Lambda +  Wl .* (Wdiag.') - Wdiag .* Wl; % Some trick to shorten the computations
  end  
  % Anti-Hermitian symmetrize at the end
  
  Lambda = Lambda - Lambda';    
end

