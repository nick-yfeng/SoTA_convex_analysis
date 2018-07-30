function [x, f, cost] = TAS_L1(y, lam0, lam1, HTH, LPF, rho, Nit)

t = 2/rho;

cost = zeros(1, Nit);       % cost function history
y = y(:);                   % convert to column vector
N = length(y);

D = @(x) diff(x); 

x = y;


for k = 1:Nit
    
    z = x - t*HTH(x - y);
    x = soft(tvd(z,N,t*lam1),t*lam0)';
    
    cost(k) = 0.5*(y - x)'*HTH(y - x) + lam1*norm(D(x),1) + lam0*norm(x,1);
end

f = LPF(y - x);     % f : low-pass component
