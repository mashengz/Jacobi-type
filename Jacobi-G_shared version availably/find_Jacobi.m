function [Psi, c,s1 ,s2] = find_Jacobi(W, i, j)
%FIND_JACOBI Summary of this function goes here
%   Detailed explanation goes here

  L = size(W, 3);
  Z = zeros(L, 3); % Matirx of Z vectors (in rows
  Z = [    squeeze(W(j,j,:) - W(i,i,:)), ...
           squeeze(W(i,j,:) + W(j,i,:)), ...
       -1i*squeeze(W(i,j,:) - W(j,i,:))];
  Gamma = 0.5 * real(Z'*Z);

  % Algorithm 4.1
  [U,S] = eig(Gamma);
  [~,i_max] = max(diag(S));
  w = U(:,i_max);
 
  if (w(1) < 0)
    w = -w;
  end
  
  c = sqrt((w(1)+1)/2);
  s1 = -w(2) / (2*c);
  s2 = -w(3) / (2*c);
  
  % Definition (2.10)
  Psi = [c         -(s1+1i*s2); ...
         s1-1i*s2    c          ];
     
end

