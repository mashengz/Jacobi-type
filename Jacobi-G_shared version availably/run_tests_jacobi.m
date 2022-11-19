function [jacobi_labels, jacobi_cost, jacobi_grad, jacobi_times, Uhs, Ahat] = ...
  run_tests_jacobi(A, Nsweeps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  n = size(A, 1);
  norm_A2 = norm(A(:),2)^2;
  sweep_size = n*(n-1)/2;
  deltas = [0.1, 0];
%   jacobi_labels = {'Jacobi-GP'; sprintf('Jacobi-G %1.1f', deltas(1));...
%     'Jacobi-C'};
  jacobi_labels = {'Jacobi-GP'; 'Jacobi-G'; 'Jacobi-C'};
  jacobi_cost = cell(length(deltas)+1,1);
  jacobi_grad = cell(length(deltas)+1,1);
  jacobi_times = cell(length(deltas)+1,1);
  Uhs = cell(length(deltas)+1,1);

   % Run Jacobi-GP 
  [Uhat, Ahat, info] = JacobiG_2P(A, deltas(1), Nsweeps);
  jacobi_cost{1} = info.iter_progress(:,1);
  jacobi_grad{1} = info.iter_progress(:,2);
  jacobi_times{1} = info.iter_times(:);
  Uhs{1} = Uhat;
  
  for i=1:length(deltas)
    [Uhat, Ahat, info] = JacobiG_2(A, deltas(i), Nsweeps);
    jacobi_cost{i+1} = info.iter_progress(:,1);
    jacobi_grad{i+1} = info.iter_progress(:,2);
    jacobi_times{i+1} = info.iter_times(:);
    Uhs{i+1} = Uhat;
  end
  
  % Check that the gradients are calculated correctly
  for i=1:length(deltas)
    [jacobi_cost{i}(end),norm_A2-matr_sumdiag2(matr_rotate(A, Uhs{i}))]
    [jacobi_grad{i}(end),norm(Lambda_2(matr_rotate(A, Uhs{i})),'fro')]
  end
end

