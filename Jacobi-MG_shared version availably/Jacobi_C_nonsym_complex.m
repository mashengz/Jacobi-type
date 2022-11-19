function [T, U, Output] = Jacobi_C_nonsym_complex(A, U_0, eps, maxiter, gtol, ftol)
%JACOBI_C_NONSYM_COMPLEX Runs the Jacobi-C algorithm for a nonsymmetric complex tensor
%   [T, U, Output] = Jacobi_C_nonsym_complex(A)
%
%   Input arguments:
%     A - complex tensor
%
%   Output arguments:
%     T - output complex tensor
%     U - cell unitary matrix
%

d     = ndims(A);
n     = size(A);

% Check the parameters
if (~exist('eps','var'))
    eps = 0; % By default, the Jacobi-C algorithm
end
if (~exist('maxiter', 'var'))
    maxiter = 1000; 
    maxsweep = 200; % sweep
end
if (~exist('gtol', 'var'))
    gtol = 1e-15;
end
if (~exist('ftol', 'var'))
    ftol = 1e-15;
end

disp('Running Jacobi-C-complex with parameters:');
fprintf('eps = %g, maxiter = %d, gtol = %g, ftol = %g \n',...
    eps, maxiter, gtol, ftol);
norm_A2 = norm(A)^2;
fprintf('Processing tensor, d = %d, ||A||^2 = %f \n', d, norm_A2);
fprintf('k | ||A||^2 - f_k \t| ||grad|| \t| (p_k,i_k,j_k) \n');
disp('------------------------------------------------------');

% Prepare the initial iteration
tic
k = 0;
if (~exist('U_0', 'var') || isempty(Q_0))
    U_k = cell(1, d);
    for i = 1 : d
        U_k{i} = eye(n(i));
    end
    W_k = A;
else
    U_k = U_0;
    W_k = tenmat_prod(A, U_k, d);
end
iter_progress = zeros(maxiter + 1, 2);
iter_sweeps   = zeros(maxiter + 1, 1);
iter_time     = zeros(maxiter + 1, 1);
f_k           = diag2sum_complex(W_k, d);
Lambda_k      = gradient_lambda_complex(W_k, d);
norm_Lambda_k = norm_cell_Lambda(Lambda_k);
% f_incr        = 1;
fprintf('%4d | %10.4e \t| %10.4e \t|       \t|       \t|       \t| \n', ...
    k, norm_A2 - f_k, norm_Lambda_k);
iter_progress(1,1) = f_k;
iter_progress(1,2) = norm_Lambda_k;
iter_progress(1,3) = 0;
iter_progress(1,4) = 0;
iter_sweeps(1,1)   = f_k;
iter_sweeps(1,2)   = norm_Lambda_k;
iter_time(1)       = toc;
sweeps             = 1;

% Execute
while (sweeps < maxsweep) %%&& (norm_Lambda_k > gtol) % && (f_incr > ftol)
    fprintf('=>*new sweeps = %d \n', sweeps);
    for p = 1 : d
        for i = 1 : (n(p) - 1)
            for j = (i + 1) : n(p)
                % Find the optimal rotation
                [x_k, f_incr] = find_maximizer_complex(W_k, i, j, p, d);
                G_k           = givens_rotation_complex(n, i, j, p, x_k);
                % Update the matrix, the tensor, the cost function value and the gradient
                U_k{p}        = U_k{p}*G_k;
                W_k           = ttm(W_k, G_k', p);
                f_k_old       = f_k;
                f_k           = f_k_old + f_incr;
                % To prevent redundant computations, the code above is
                % equivalent to, but lower computations than, the code below
                %                     f_k           = diag2sum_complex(W_k, d);
                Lambda_k      = gradient_lambda_complex(W_k, d);
                norm_Lambda_k = norm_cell_Lambda(Lambda_k);
                k             = k + 1;
                fprintf('%4d | %10.4e \t| %10.4e \t| (%2d,%2d,%2d) \n', ...
                    k, norm_A2 - f_k, norm_Lambda_k, p, i, j);
                iter_progress(k+1, 1) = f_k;
                iter_progress(k+1, 2) = norm_Lambda_k;
                iter_progress(k+1, 3) = x_k(1);
                iter_progress(k+1, 4) = f_incr;
                iter_time(k+1)        = toc;
            end
        end
    end
    sweeps = sweeps + 1;
    iter_sweeps(sweeps+1,1) = f_k;
    iter_sweeps(sweeps+1,2) = norm_Lambda_k;
%     if (k >= maxiter)
%         break;
%     end
end
U                    = U_k;
T                    = W_k;
Output.iter          = k;
Output.sweeps        = sweeps;
Output.iter_sweeps   = iter_sweeps;
Output.iter_time     = iter_time;
Output.iter_progress = iter_progress;
if (sweeps >= maxiter)
    fprintf('Exceed the maximum iteration number: sweeps is greater than %d \n', ...
        maxiter);
else
    fprintf('Gradient norm tolerance reached: gtol = %g \n', ...
        gtol);
end
fprintf('Total time is %f [s] \n', ...
    Output.iter_time(k));
end

