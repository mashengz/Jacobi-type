function [Gamma] = Gamma_2(W)
%GAMMA_2 Summary of this function goes here
%   Returns the Gamma matrix for joint diagonalization of L 2x2 matrices

L = size(W, 3);
Z = zeros(L, 3); % Matirx of Z vectors (in rows

Z = [    squeeze(W(2,2,:) - W(1,1,:)), ...
         squeeze(W(1,2,:) + W(2,1,:)), ...
     -1i*squeeze(W(1,2,:) - W(2,1,:))];
Gamma = 0.5 * (norm(squeeze(W(1,1,:)+ W(2,2,:)),2) * eye(3) + real(Z'*Z));
end