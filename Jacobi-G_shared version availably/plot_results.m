function [f_cost,f_grad] = plot_results(labels, times, cost, grad, sty)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p = 'symm';
f_cost = figure
for i=1:length(labels)
  semilogy(times{i}-times{i}(1), cost{i}, sty{i}); hold on;
end  
xlabel('Time (s)');
ylabel('$\sum_{\ell}\|\mathbf{A}^{(\ell)}\|^2 - g(\mathbf{U}_k)$','Interpreter','Latex','FontSize', 10);
% title('Cost function value'); 
legend(labels{:});
xlim([0,0.15]);ylim([10^(-0.5),10^(2.5)]);
figure_1 = ['ex_' num2str(p) '_fun_time'  '.eps' ];
saveas(gcf, figure_1, 'epsc');
hold off;

f_grad = figure
for i=1:length(labels)
  semilogy(times{i}-times{i}(1), grad{i}, sty{i}); hold on;
end
% title('Norm of the gradient'); 
xlabel('Time (s)');
ylabel('Norm of the gradient','Interpreter','Latex','FontSize', 10);
legend(labels{:});
xlim([0,0.15]);ylim([10^(-20),10^(5)]);
figure_2 = ['ex_' num2str(p) '_grad_time'  '.eps' ];
saveas(gcf, figure_2, 'epsc');
hold off;
end

