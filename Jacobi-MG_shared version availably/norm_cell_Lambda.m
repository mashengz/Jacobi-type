function [g] = norm_cell_Lambda(A)
%NORM_CELL_LAMBDA Returns the Frobenius norm of the cell matrix
%   [g] = norm_cell_Lambda(A)
%
%   Input argument:
%     A - cell matrix, input argument of type 'cell'
%
%   Output argument:
%     g - the Frobenius norm of the cell matrix
%

m = size(A, 2);
g = 0;
for i = 1 : m
    g = g + frob(A{i})^2;
end
g = sqrt(g);
end
    