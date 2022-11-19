function [info_sd, info_cg, info_bb, T, X, Y, Z] = manopt_nonsym_3rd_complex_r(A,r,options)
%MANOPT_NONSYM_3rd_COMPLEX Runs the Riemannian algorithm for a 3-rd nonsymmetric complex tensor
%   [info_sd, info_cg, info_bb, T, X, Y, Z] = manopt_nonsym_3rd_complex_r(A)
%
%   Input arguments:
%     A - complex tensor
%     r - specified low dimensions
%
%   Output arguments:
%     info - output parameters for different Riemannian algorithms
%     T - output complex tensor
%     X, Y, Z - output approximate unitary matrices/optimal solutions
%

% Verify that Manopt was indeed added to the Matlab path
if isempty(which('stiefelcomplexfactory'))
    error(['You should first add Manopt to the Matlab path.\n' ...
        'Please run importmanopt.']);
end

% Generate the problem data
d     = ndims(A);
n     = size(A);
r     = r*ones(1,d);
n_min = min(n);
MM    = cell(1, d);
for i = 1 : d
    MM{i} = [eye(r(i)),zeros(r(i),n(i)-r(i));...
        zeros(n(i)-r(i),r(i)),zeros(n(i)-r(i),n(i)-r(i))];
end
MX = MM{1};
MY = MM{2};
MZ = MM{3};

% Create the problem structure
tuple.X = stiefelcomplexfactory(n(1), n(1));
tuple.Y = stiefelcomplexfactory(n(2), n(2));
tuple.Z = stiefelcomplexfactory(n(3), n(3));
M = productmanifold(tuple);

% Define the problem cost function and its gradient
problem.M = M;
problem.cost  = @cost_function;
problem.egrad = @egrad_function;

% Create the problem cost function
    function [f] = cost_function(y)
        X = y.X;
        Y = y.Y;
        Z = y.Z;
        W = ttm(A, {MX*X',MY*Y',MZ*Z'}, [1 2 3]);
        f = 0;
        for k = 1 : r
            f = f + abs(W(k,k,k))^2;
        end
        norm_A2 = frob(A.data)^2;
        f       = norm_A2 - f;
    end

% Create the probelm Euclidean gradient function
    function [g] = egrad_function(y)
        X = y.X;
        Y = y.Y;
        Z = y.Z;
        W = ttm(A, {MX*X',MY*Y',MZ*Z'}, [1 2 3]);
        v = zeros(n_min, 1);
        matr_arg      = cell(1, d);
        [matr_arg{:}] = deal(':');
        for l = 1 : n_min
            [matr_arg{1:d}] = deal(l);
            v(l) = W(matr_arg{:});
        end
        D = tendiag(v,n);
        W1 = tens2mat(D.data, 1);
        W2 = tens2mat(D.data, 2);
        W3 = tens2mat(D.data, 3);
        BX = ttm(A, {MY*Y',MZ*Z'}, [2 3]);
        B1 = tens2mat(BX.data, 1);
        lambdaX = (B1*W1')*MX;
        g.X = - 2*lambdaX;
        BY = ttm(A, {MX*X',MZ*Z'}, [1 3]);
        B2 = tens2mat(BY.data, 2);
        lambdaY = (B2*W2')*MY;
        g.Y = - 2*lambdaY;
        BZ = ttm(A, {MX*X',MY*Y'}, [1 2]);
        B3 = tens2mat(BZ.data, 3);
        lambdaZ = (B3*W3')*MZ;
        g.Z = - 2*lambdaZ;
    end


% Generate an initial guess
% Random initial guess
% x0 = problem.M.rand();
% Identity matrix as an initial guess
E = cell(1, d);
for j = 1 : d
    E{j} = eye(n(j));
end
scores = {'X';'Y';'Z'};
x0 = cell2struct(E, scores, 2);

% Solve
disp('Running Steepest descent algorithm');
[x_sd, xcost_sd, info_sd] = steepestdescent(problem, x0, options);
disp('------------------------------------------------------------------');
disp('Running Conjugate gradient algorithm');
[x_cg, xcost_cg, info_cg] = conjugategradient(problem, x0, options);
disp('------------------------------------------------------------------');
disp('Running Barzilai Borwein algorithm');
[x_bb, xcost_bb, info_bb] = barzilaiborwein(problem, x0, options);
disp('------------------------------------------------------------------');
T.sd = ttm(A, {x_sd.X', x_sd.Y', x_sd.Z'}, [1 2 3]);
T.cg = ttm(A, {x_cg.X', x_cg.Y', x_cg.Z'}, [1 2 3]);
T.bb = ttm(A, {x_bb.X', x_bb.Y', x_bb.Z'}, [1 2 3]);
end
