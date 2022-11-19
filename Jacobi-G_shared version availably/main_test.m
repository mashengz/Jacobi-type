%% Main test reproducing example from the article
clear;
clc;
close all;

%% Diagonal example
% How it was generated
n = 10; L= 8;
% D = zeros(n,n,L);
% 
% In = eye(n);
% for i=1:L
%   D(:,:,i) = diag(2*randn(n,1)); 
% end
% X = randn(n);
% [Ust,~,~] = svd(X);
% E = tensor(2*randn(n,n,L));
% A = matr_rotate(D, Ust) + E.data/norm(E);
% save('diag2_mat.mat','A');
load('diag2_mat.mat','A');  % diag2_mat.mat is an example in SLN, Jacobi_Poly. Final
% 
norm_A2 = 0;
for i=1:L
    norm_A2 = norm_A2 + norm(tensor(A(:,:,i)))^2;
end
diag_squares_sum_A = 0;
for i=1:L
  diag_squares_sum_A = diag_squares_sum_A + norm(diag(A(:,:,i)))^2;
end
percentage_A = diag_squares_sum_A/norm_A2;

%% Optimize with various algorithms
Nsweeps = 20;
[jacobi_labels, jacobi_cost, jacobi_grad, jacobi_times, Uhs, Ahat] = ...
    run_tests_jacobi(A, Nsweeps);
[manopt_solvers,manopt_labels,manopt_cost, ...
    manopt_grad,manopt_info,manopt_times,manopt_Uh] = run_tests_manopt(A);

%% Compute percentage
percentage_Ahat = Ahat/norm_A2;

%% Plot results
[f3,f4] = plot_results([ manopt_labels;jacobi_labels],[manopt_times;jacobi_times], ...
             [manopt_cost;jacobi_cost],[manopt_grad;jacobi_grad], ...
             {'-m'; '-g'; '-r';'-b'; '-k'; '-c'});





    