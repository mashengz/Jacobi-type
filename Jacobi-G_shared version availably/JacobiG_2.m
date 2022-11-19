function [Uhat, Ahat, info] = JacobiG_2(A, delta_rel, maxsweep, U0, grad_eps, x_eps)
%JACOBI_G_2 Run the Jacobi-G algorithm for nxnxL  
% if delta = 0, then it is Jacobi-C
  n = size(A,1);
  d = 2;
  L = size(A,3);

  % Check the parameters
  if (~exist('delta_rel','var')), delta_rel = 0; end 
  % By default, the Jacobi-C algorithm.   
  if (~exist('maxsweep', 'var')), maxsweep = 100; end
  if (~exist('grad_threshold')), grad_eps = 1e-15; end
  if (~exist('x_threshold')), x_eps = 1e-16; end
  if (~exist('opt_disp')), opt_disp = 0; end

  delta = delta_rel * sqrt(2) / n;
  
  disp('Running Jacobi-G with parameters:');
  fprintf(['delta = %g, maxsweep = %d, grad_eps = %g, x_eps = %g\n'],...
       delta, maxsweep, grad_eps, x_eps);
  norm_A2 = norm(A(:))^2;
  fprintf('Processing tensor, d = %d, n = %d, ||A||^2 = %f\n', ...
               d, n, norm(A(:))^2);  
  if (opt_disp ~= 0)
    fprintf(['   k | f(U_k) - ||A||^2\t| ||ProjGrad||  \t|' ...
                  '(f_k-f_{k-1})/(h''(0)|x_k|)\n']);     
    disp(['--------------------------------------------------------------'...
         '----------------------------']);
  end
  sweep_size = (n * (n-1)) /2;
   
  tic
  
  % Prepare the zero-th iteration
  W_k = zeros(n,n,L);  
  if (~exist('U0', 'var') || isempty(U0))
    U_k = eye(n); W_k = A;
  else  
    U_k = U0; W_k = matr_rotate(A, U0);
  end
  
  
  iter_pairs = zeros(maxsweep+1, 2);
  iter_progress = zeros(maxsweep+1, 2);
  iter_times =  zeros(maxsweep+1,1);
  
  k = 0;
  f_k = matr_sumdiag2(W_k);
  Lambda_k = Lambda_2(W_k);
  norm_Lambda_k = norm(Lambda_k, 'fro');
  if (opt_disp ~= 0)
    fprintf('%4d | %10.5e \t| %10.5e \t|         \t|\n', ...
       k, norm_A2-f_k, norm_Lambda_k);
  end   
  iter_pairs(1,:) = [0,0];
  iter_progress(1,:) = [norm_A2- f_k,norm_Lambda_k];
  iter_times(1) = toc;
  
  iters = 0;
  for k=1:maxsweep
    for i=1:n
      for j=i+1:n
        h0 = sqrt(2)* abs(Lambda_k(i,j));
        if (abs(h0) >= delta * norm_Lambda_k) % Check the inequality
          iters = iters + 1; 
          [Psi_k,c_k,~,~] = find_Jacobi(W_k, i, j);
          G_k = givens_complex(n, i, j, Psi_k);
          W_k = matr_rotate(W_k, G_k); 
          U_k(:, [i,j]) =  U_k(:, [i,j]) * Psi_k;
       
          Lambda_k = Lambda_2(W_k);
          if (delta  > 0)
            norm_Lambda_k = norm(Lambda_k, 'fro');
          end
        end
      end
    end  
    if (delta == 0)
      norm_Lambda_k = norm(Lambda_k, 'fro');
    end
    f_k1 = f_k;
    f_k = matr_sumdiag2(W_k);
   
    if (opt_disp ~= 0)
      fprintf('%4d | %10.5e \t| %10.5e  \t|%10.5e \t|\n', ...
        k, norm_A2- f_k, norm_Lambda_k, (f_k-f_k1)/abs(h0*c_k));
    end
    iter_pairs(k+1,:) = [i,j];
    iter_progress(k+1,:) = [norm_A2-f_k , norm_Lambda_k];
    iter_times(k+1) = toc;
  end
  
  info.iter = iters;
  info.sweeps = k;
  info.iter_pairs = iter_pairs;
  info.iter_progress = iter_progress;
  info.iter_times = iter_times;
  Uhat = U_k;
  Ahat = f_k;
end

