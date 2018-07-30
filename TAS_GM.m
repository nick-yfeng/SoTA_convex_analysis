function [x, u, f, cost] = TAS_GM(y, lam0, lam1, HTH, LPF, gamma, rho, Nit, Innit)


t = 2/rho;

cost = zeros(1, Nit);       % cost function history
y = y(:);                   % convert to column vector
N = length(y);

D = @(x) diff(x); 


Innercost = zeros(1, Innit);

x = y;
u_prev = x;

for k = 1:Nit
    
    u = u_prev;
    for i = 1:Innit
        
        w1 = u - t*HTH(u - x);
        u = soft(tvd(w1,N,t*lam1/gamma),t*lam0/gamma)';

        Innercost(i) = 0.5*gamma*(u-x)'*HTH(u-x) + lam1*norm(D(u),1) + lam0*norm(u,1);  
    
    end
    u_prev = u;
%     figure(1)
%     plot(Innercost)
%     drawnow
    
    
    w2 = x - t*HTH((1-gamma)*x + gamma*u - y);
    x = soft(tvd(w2,N,t*lam1),t*lam0)';
    
    
    cost(k) = 0.5*(y-x)'*HTH(y-x) - 0.5*gamma*(u-x)'*HTH(u-x) + lam1*(norm(D(x),1) - norm(D(u),1)) + lam0*(norm(x,1) - norm(u,1));
end

                    
f = LPF(y - x);     % f : low-pass component
