function [info_sd, info_cg, info_bb, T, X, Y, Z, Q] = manopt_nonsym_4th_complex(A,options)
%MANOPT_NONSYM_4th_COMPLEX Runs the Riemannian algorithm for a 4-th nonsymmetric complex tensor
%   [info_sd, info_cg, info_bb, T, X, Y, Z, Q] = manopt_nonsym_4th_complex(A)
%
%   Input arguments:
%     A - complex tensor
%
%   Output arguments:
%     info - output parameters for different Riemannian algorithms
%     T - output complex tensor
%     X, Y, Z, Q - output approximate unitary matrices/optimal solutions
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
tuple.Q = stiefelcomplexfactory(n(4), n(4));
M = productmanifold(tuple);

% Define the problem cost function and its gradient
problem.M = M;
problem.cost  = @cost_function;
problem.egrad = @egrad_function;

% Create the probelm cost function
    function [f] = cost_function(y)
        X = y.X;
        Y = y.Y;
        Z = y.Z;
        Q = y.Q;
        W = ttm(A, {X',Y',Z',Q'}, [1 2 3 4]);
        f = 0;
        for i = 1 : n_min
            f = f + abs(W(i,i,i,i))^2;
        end
        norm_A2 = norm(A)^2;
        f       = norm_A2 - f;
    end

% Create the probelm Euclidean gradient function
    function [g] = egrad_function(y)
        X = y.X;
        Y = y.Y;
        Z = y.Z;
        Q = y.Q;
        W = ttm(A, {X',Y',Z',Q'}, [1 2 3 4]);
        lambdaX = zeros(n(1), n(1));
        for l = 1 : n(1)
            for k = 1 : n(1)
                if (l <= n_min)
                    lambdaX(l,l) = abs(W([l,l,l,l]))^2;
                else
                    lambdaX(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaX(l,k) = W([l,k,k,k])*conj(W([k,k,k,k]));
                end
                if (l <= n_min)
                    lambdaX(k,l) = W([k,l,l,l])*conj(W([l,l,l,l]));
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
                    lambdaY(l,l) = abs(W([l,l,l,l]))^2;
                else
                    lambdaY(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaY(l,k) = W([k,l,k,k])*conj(W([k,k,k,k]));
                end
                if (l <= n_min)
                    lambdaY(k,l) = W([l,k,l,l])*conj(W([l,l,l,l]));
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
                    lambdaZ(l,l) = abs(W([l,l,l,l]))^2;
                else
                    lambdaZ(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaZ(l,k) = W([k,k,l,k])*conj(W([k,k,k,k]));
                end
                if (l <= n_min)
                    lambdaZ(k,l) = W([l,l,k,l])*conj(W([l,l,l,l]));
                end
                if (l > n_min) && (k > n_min)
                    lambdaZ(l,k) = 0;
                    lambdaZ(k,l) = 0;
                end
            end
        end
        g.Z = - 2*Z*lambdaZ;
        lambdaQ = zeros(n(4), n(4));
        for l = 1 : n(4)
            for k = 1 : n(4)
                if (l <= n_min)
                    lambdaQ(l,l) = abs(W([l,l,l,l]))^2;
                else
                    lambdaQ(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaQ(l,k) = W([k,k,k,l])*conj(W([k,k,k,k]));
                end
                if (l <= n_min)
                    lambdaQ(k,l) = W([l,l,l,k])*conj(W([l,l,l,l]));
                end
                if (l > n_min) && (k > n_min)
                    lambdaQ(l,k) = 0;
                    lambdaQ(k,l) = 0;
                end
            end
        end
        g.Q = - 2*Q*lambdaQ;
    end

% Generate an initial guess
% Random initial guess
% x0 = problem.M.rand();
% Identity matrix as an initial guess
E = cell(1, d);
for j = 1 : d
    E{j} = eye(n(j));
end
scores = {'X';'Y';'Z';'Q'};
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
T.sd = ttm(A, {x_sd.X', x_sd.Y', x_sd.Z', x_sd.Q'}, [1 2 3 4]);
T.cg = ttm(A, {x_cg.X', x_cg.Y', x_cg.Z', x_cg.Q'}, [1 2 3 4]);
T.bb = ttm(A, {x_bb.X', x_bb.Y', x_bb.Z', x_bb.Q'}, [1 2 3 4]);
end
