function [info_sd, info_cg, info_bb, T, X, Y] = manopt_nonsym_2nd_complex_r(A,r,options)
%MANOPT_NONSYM_2nd_COMPLEX Runs the Riemannian algorithm for a 2-nd nonsymmetric complex tensor
%   [info_sd, info_cg, info_bb, T, X, Y] = manopt_nonsym_2nd_complex_r(A)
%
%   Input arguments:
%     A - complex tensor
%     r - specified low dimensions
%
%   Output arguments:
%     info - output parameters for different Riemannian algorithms
%     T - output complex tensor
%     X, Y - output approximate unitary matrices/optimal solutions
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

% Create the problem structure
tuple.X = stiefelcomplexfactory(n(1), n(1));
tuple.Y = stiefelcomplexfactory(n(2), n(2));
M = productmanifold(tuple);

% Define the problem cost function and its gradient
problem.M = M;
problem.cost  = @cost_function;
problem.egrad = @egrad_function;

% Create the problem cost function
    function [f] = cost_function(y)
        X = y.X;
        Y = y.Y;
        W = ttm(A, {X',Y'}, [1 2]);
        f = 0;
        for i = 1 : r
            f = f + abs(W(i,i))^2;
        end
        norm_A2 = norm(A)^2;
        f       = norm_A2 - f;
    end

% Create the probelm Euclidean gradient function
    function [g] = egrad_function(y)
        X = y.X;
        Y = y.Y;
        W = ttm(A, {X',Y'}, [1 2]);
        lambdaX = zeros(n(1), n(1));
        for l = 1 : n(1)
            for k = 1 : n(1)
                if (l <= n_min)
                    lambdaX(l,l) = abs(W([l,l]))^2;
                else
                    lambdaX(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaX(l,k) = W([l,k])*conj((W([k,k]));
                end
                if (l <= n_min)
                    lambdaX(k,l) = W([k,l])*conj(W([l,l]));
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
                    lambdaY(l,l) = abs(W([l,l]))^2;
                else
                    lambdaY(l,l) = 0;
                end
                if (k <= n_min) && (l < k) && (l <= n_min)
                    lambdaY(l,k) = W([k,l])*conj(W([k,k]));
                end
                if (l <= n_min)
                    lambdaY(k,l) = W([l,k])*conj(W([l,l]));
                end
                if (l > n_min) && (k > n_min)
                    lambdaY(l,k) = 0;
                    lambdaY(k,l) = 0;
                end
            end
        end
        g.Y = - 2*Y*lambdaY;
    end


% Generate an initial guess
% Random initial guess
% x0 = problem.M.rand();
% Identity matrix as an initial guess
E = cell(1, d);
for j = 1 : d
    E{j} = eye(n(j));
end
scores = {'X';'Y'};
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
T.sd = ttm(A, {x_sd.X', x_sd.Y'}, [1 2]);
T.cg = ttm(A, {x_cg.X', x_cg.Y'}, [1 2]);
T.bb = ttm(A, {x_bb.X', x_bb.Y'}, [1 2]);
end
