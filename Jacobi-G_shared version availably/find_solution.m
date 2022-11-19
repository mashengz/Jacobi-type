function [x] = find_solution(Gamma,eps,e)

[U,S] = eig(Gamma);
[~,i_max] = max(diag(S));
x0 = U(:,i_max);
x0 = x0/norm(x0);
lambda0 = -2*S(i_max,i_max);
epsilon = 10^(-5);
N = 3000;
for k = 1:N 
    A = [Gamma+lambda0*eye(3,3), x0; -2*x0', 0];
    F = [(Gamma+lambda0*eye(3,3))*x0-eps*e;1-x0'*x0];
    d = A\(-F);
    x = x0+d(1:3);
    lambda0 = lambda0+d(4);
    norm_d = norm(F);
    if (norm_d < epsilon)
        return
    else
        x0 = x;
        if  k == N
            warning('error')
        end
    end
end
        
    
