function [x,u] = LF(H,K,POS,PLOT)
h = 2^-H;
k = 2^-K;
lambda = k/h;

x = -1.25:h:1.25;
xl = length(x);
u = zeros(xl,1);
if POS == 1
    u(x >= 0) = 1;
else
    u(x < 0) = 1;
end

Dt = sparse(diag(ones(xl-1,1),1));
Db = sparse(diag(ones(xl-1,1),-1));
Dt([1,end],:) = 0;
Db([1,end],:) = 0;
A1 = Dt + Db;
A2 = Dt - Db;
A1(1,1) = 1;
A1(end,end) = 1;

for i = 1:(1/k)
    u = A1*u - lambda*(A2*(1/2*u.^2));
    if PLOT
        plot(x,u)
        axis([min(x) max(x) 0 1])
        drawnow
    end
end
end