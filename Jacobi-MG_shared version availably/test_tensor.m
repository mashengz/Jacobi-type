function [A, r] = test_tensor(ind)
%TEST_TENSOR Creates a tensor
%   [A, r] = test_tensor(ind)
%
%   Input argument:
%     ind - index
%
%   Output argument:
%     A - test tensor
%     r - specified low dimensions (r == 0 means this step can be ignored)
%
% Written by sz

%-------------------------- Real tensor -----------------------------------
if (ind == 1) 
    A(:,:,1) = [8 8 3;10 5 7;10 5 4];
    A(:,:,2) = [10 8 10;8 3 7;5 5 3];
    A(:,:,3) = [9 3 4;7 7 6;2 7 5];
    A = tensor(A);
	r = 0;
end
%-------------------------- Complex tensor --------------------------------
if (ind == 1i)  
%     n = [5 5 6 6];
    Construct_tensor = load('ex_0+1i.mat');
    A = Construct_tensor.A;
	r = 0;
end

if (ind == 2i)  
%     n = [3 3 3 3];
    Construct_tensor = load('ex_0+2i.mat');
    A = Construct_tensor.A;
	r = 2;
end
end

