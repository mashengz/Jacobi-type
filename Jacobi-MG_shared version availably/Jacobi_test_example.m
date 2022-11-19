% function Jacobi_test_example()
%
% Written by sz

clear;
clc;
close all;

% Generate the problem data: tensor
ind     = 1i;
[A, r]  = test_tensor(ind);
d       = ndims(A);
norm_A2 = norm(A)^2;

% Prepare the functions
if (r == 0)
    JacobiG_nonsym = @JacobiG_nonsym_complex;
    JacobiC_nonsym = @Jacobi_C_nonsym_complex;
    if (d == 2)
        manopt_nonsym = @manopt_nonsym_2nd_complex;
    elseif (d == 3)
        manopt_nonsym = @manopt_nonsym_3rd_complex;
    elseif (d == 4)
        manopt_nonsym = @manopt_nonsym_4th_complex;
    else
        error('Unsupported order of the tensor (only 2, 3, 4 are supported)');
    end
else
    JacobiG_nonsym = @JacobiG_nonsym_complex_r;
    JacobiC_nonsym = @Jacobi_C_nonsym_complex_r;
    if (d == 2)
        manopt_nonsym = @manopt_nonsym_2nd_complex_r;
    elseif (d == 3)
        manopt_nonsym = @manopt_nonsym_3rd_complex_r;
    elseif (d == 4)
        manopt_nonsym = @manopt_nonsym_4th_complex_r;
    else
        error('Unsupported order of the tensor (only 2, 3, 4 are supported)');
    end
end

% Calculate the percentage of both the sums of squares of the diagonals and the Fro-norm of tensor A
diag_squares_sum_A = 0;
for i = 1 : min(size(A))
    diag_squares_sum_A = diag_squares_sum_A + abs(A(i.*ones(1,d)))^2;
end
percentage_A = diag_squares_sum_A/norm_A2;

% Execute some Riemannian algorithms and Jacobi-G nonsymmetric algorithm
options.tolgradnorm = 1e-15;
options.minstepsize = 0;
% options.maxiter = 1000;
if (r == 0)
    [info_sd, info_cg, info_bb, T_r] = manopt_nonsym(A,options);
    [T, U, Output]                   = JacobiG_nonsym(A);
    [T_c, U_c, Output_c]             = JacobiC_nonsym(A);
else
    [info_sd, info_cg, info_bb, T_r] = manopt_nonsym(A,r,options);
    [T, U, Output]                   = JacobiG_nonsym(A,r);
    [T_c, U_c, Output_c]             = JacobiC_nonsym(A,r);
end

% Calculate the percentage of both the sums of squares of the diagonals and the Fro-norm of tensor T
diag_squares_sum    = 0;
diag_squares_sum_sd = 0;
diag_squares_sum_cg = 0;
diag_squares_sum_bb = 0;
diag_squares_sum_c  = 0;
if (r == 0)
    r = min(size(T));
end
for i = 1 : r
    diag_squares_sum    = diag_squares_sum + abs(T(i.*ones(1,d)))^2;
    diag_squares_sum_sd = diag_squares_sum_sd + abs(T_r.sd(i.*ones(1,d)))^2;
    diag_squares_sum_cg = diag_squares_sum_cg + abs(T_r.cg(i.*ones(1,d)))^2;
    diag_squares_sum_bb = diag_squares_sum_bb + abs(T_r.bb(i.*ones(1,d)))^2;
    diag_squares_sum_c  = diag_squares_sum_c + abs(T_c(i.*ones(1,d)))^2;
end
percentage    = diag_squares_sum/norm_A2;
percentage_sd = diag_squares_sum_sd/norm_A2;
percentage_cg = diag_squares_sum_cg/norm_A2;
percentage_bb = diag_squares_sum_bb/norm_A2;
percentage_c  = diag_squares_sum_c/norm_A2;

% Display some figures
iter    = Output.sweeps;
iter_c  = Output_c.sweeps;
time    = Output.iter;
time_c  = Output_c.iter;
time_sd = [info_sd.time];
time_cg = [info_cg.time];
time_bb = [info_bb.time];

figure;
semilogy(time_sd - time_sd(1), [info_sd.cost], '-m');
hold on;
semilogy(time_cg - time_cg(1), [info_cg.cost], '-g');
hold on;
semilogy(time_bb - time_bb(1), [info_bb.cost], '-r');
hold on;
semilogy(Output.iter_time(1:time) - Output.iter_time(1), norm_A2*ones(time, 1) - Output.iter_progress(1:time, 1), '-b');
hold on;
semilogy(Output_c.iter_time(1:time_c) - Output_c.iter_time(1), norm_A2*ones(time_c, 1) - Output_c.iter_progress(1:time_c, 1), '-k');
xlabel('Time (s)');
ylabel('$\|\mathcal{A}\|^2 - f(\mathbf{\omega}_k)$','Interpreter','Latex','FontSize', 10);
legend('Steepest descent','Conjugate gradient',...
    'Barzilai Borwein','Jacobi-MG','Jacobi-MC')
% xlim([0, 1]); % ind = 1
if (ind == 1i)
    xlim([0, 6]); % ind = 1i
else
    xlim([0, 2]);  % ind = 2i
    ylim([10^(1.2),1.2*max(norm_A2*ones(time, 1) - Output.iter_progress(1:time, 1))]);  % ind = 2i
end
figure_1 = ['ex_' num2str(ind) '_fun_time'  '.eps' ];
saveas(gcf, figure_1, 'epsc');


figure;
semilogy(time_sd - time_sd(1), [info_sd.gradnorm], '-m');
hold on;
semilogy(time_cg - time_cg(1), [info_cg.gradnorm], '-g');
hold on;
semilogy(time_bb - time_bb(1), [info_bb.gradnorm], '-r');
hold on;
semilogy(Output.iter_time(1:time) - Output.iter_time(1), Output.iter_progress(1:time, 2), '-b');
hold on;
semilogy(Output_c.iter_time(1:time_c) - Output_c.iter_time(1), Output_c.iter_progress(1:time_c, 2), '-k');
xlabel('Time (s)');
ylabel('Norm of the gradient','Interpreter','Latex','FontSize', 10);
legend('Steepest descent','Conjugate gradient',...
    'Barzilai Borwein','Jacobi-MG','Jacobi-MC')
% xlim([0, 4]);  % ind = 1
if (ind == 1i)
    xlim([0, 18]); % ind = 1i
    ylim([1e-17, 10*max(Output.iter_progress(1:time, 2))]);  % ind = 1i
else
    xlim([0, 8]); % ind = 2i
end
figure_2 = ['ex_' num2str(ind) '_grad_time'  '.eps' ];
saveas(gcf, figure_2, 'epsc');

