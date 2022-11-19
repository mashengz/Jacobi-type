function [manopt_solvers,manopt_labels,manopt_cost, ...
manopt_grad,manopt_info,manopt_times,manopt_Uh] = run_tests_manopt(A)
  n = size(A,1);
%   norm_A2 = norm(A(:),2)^2;
  % Create the problem structure.
  manifold = stiefelcomplexfactory(n,n);
  problem.M = manifold;
  
  norm_A2 = norm(A(:))^2;
  
  % Define the problem cost function and its gradient.
  problem.cost = @(U) norm_A2-matr_sumdiag2(matr_rotate(A, U));
  problem.grad = @(U) -(U * Lambda_2(matr_rotate(A, U)));
  
  % Numerically check gradient and Hessian consistency.
  figure;
  checkgradient(problem);
  
  %opt.max_iter = Nsweeps;
%   opt = struct ('tolgradnorm', 1e-15, 'minstepsize', 1e-14, 'verbosity', 1);
  opt = struct ('tolgradnorm', 0, 'minstepsize', 0, 'verbosity', 1, 'maxiter', 1000);
  manopt_solvers = {'steepestdescent'; 'conjugategradient'; 'barzilaiborwein'};
  manopt_labels = {'Steepest descent'; 'Conjugate gradient'; 'Barzilai Borwein'};
  manopt_cost = cell(length(manopt_solvers),1);
  manopt_grad = cell(length(manopt_solvers),1);
  manopt_info = cell(length(manopt_solvers),1);
  manopt_times = cell(length(manopt_solvers),1);
  manopt_Uh = cell(length(manopt_solvers),1);
  
  
  for i=1:size(manopt_solvers,1)
    eval(['[manopt_Uh{i}, manopt_cost{i}, manopt_info{i}, ~] = '...
      manopt_solvers{i} '(problem, eye(n),opt);'])
    manopt_cost{i} = [manopt_info{i}.cost]';
    manopt_grad{i} = [manopt_info{i}.gradnorm]';
    manopt_times{i} = [manopt_info{i}.time]';
    % Check that the gradients are calculated correctly
    [manopt_grad{i}(end),norm(Lambda_2(matr_rotate(A, manopt_Uh{i})),'fro')]
  end
end

