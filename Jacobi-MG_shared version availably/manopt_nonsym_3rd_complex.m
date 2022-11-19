function [info_sd, info_cg, info_bb, T, X, Y, Z] = manopt_nonsym_3rd_complex(A,options)
%MANOPT_NONSYM_3rd_COMPLEX Runs the Riemannian algorithm for a 3-rd nonsymmetric complex tensor
%   [info_sd, info_cg, info_bb, T, X, Y, Z] = manopt_nonsym_3rd_complex(A)
%
%   Input arguments:
%     A - complex tensor
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
d = ndims(A);
n = size(A);
n_min = min(n);

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
        W = ttm(A, {X',Y',Z'}, [1 2 3]);
        f = 0;
        for i = 1 : n_min
            f = f + abs(W(i,i,i))^2;
        end
        norm_A2 = norm(A)^2;
        f       = norm_A2 - f;
    end

% Create the probelm Euclidean gradient function
    function [g] = egrad_function(y)
        X = y.X;
        Y = y.Y;
        Z = y.Z;
        W = ttm(A, {X',Y',Z'}, [1 2 3]);
        lambdaX = zeros(n(1), n(1));
        for l = 1 : n(1)
            for k = 1 : n(1)
                if (l <= n_min)
                    lambdaX(l,l) = abs(W([l,l,l]))^2;
                else
                    lambdaX(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaX(l,k) = W([l,k,k])*conj(W([k,k,k]));
                end
                if (l <= n_min)
                    lambdaX(k,l) = W([k,l,l])*conj(W([l,l,l]));
                end
                if (l > n_min) && (k > n_min)
                    lambdaX(l,k) = 0;
                    lambdaX(k,l) = 0;
                end
            end
        end
        g.X = - 2*X*lambdaX;
        lambdaY = zeros(n(2), n(2));
        for l = 1 : n(2)
            for k = 1 : n(2)
                if (l <= n_min)
                    lambdaY(l,l) = abs(W([l,l,l]))^2;
                else
                    lambdaY(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaY(l,k) = W([k,l,k])*conj(W([k,k,k]));
                end
                if (l <= n_min)
                    lambdaY(k,l) = W([l,k,l])*conj(W([l,l,l]));
                end
                if (l > n_min) && (k > n_min)
                    lambdaY(l,k) = 0;
                    lambdaY(k,l) = 0;
                end
            end
        end
        g.Y = - 2*Y*lambdaY;
        lambdaZ = zeros(n(3), n(3));
        for l = 1 : n(3)
            for k = 1 : n(3)
                if (l <= n_min)
                    lambdaZ(l,l) = abs(W([l,l,l]))^2;
                else
                    lambdaZ(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaZ(l,k) = W([k,k,l])*conj(W([k,k,k]));
                end
                if (l <= n_min)
                    lambdaZ(k,l) = W([l,l,k])*conj(W([l,l,l]));
                end
                if (l > n_min) && (k > n_min)
                    lambdaZ(l,k) = 0;
                    lambdaZ(k,l) = 0;
                end
            end
        end
        g.Z = - 2*Z*lambdaZ;
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
