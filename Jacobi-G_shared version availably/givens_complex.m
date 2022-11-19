function [G] = givens_complex(n, i, j, Psi)
%GIVENS_ROTATION Creates a  complex Givens rotation based on Psi
%   [G] = givens_rotation(n, i, j, Psi)
%
%   Input arguments:
%     n - size of the matrix
%     i,j - pair of indices
%     Psi - matrix
%
%   Output arguments:
%     G - the givens matrix
  G = speye(n);
  G([i j],[i,j]) = Psi;
end

