function [x,u] = LW(H,K,POS,e,PLOT)
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
Dc = sparse(diag(ones(xl,1),0));
Dt([1,end],:) = 0;
Db([1,end],:) = 0;
Dc([1,end],:) = 0;
Ac = Dt/2 - Db/2;
Af = Dt/2 - Dc/2;
Ab = Dc/2 - Db/2;
Avf = Dt/2 + Dc/2;
Avb = Dc/2 + Db/2;
A2c = Dt + Db - 2*Dc;
Avf(1,1) = 1;
Avb(1,1) = 1;
Avf(end,end) = 1;
Avb(end,end) = 1;

for i = 1:(1/k)
    f = 1/2*u.^2;
    al = Avb*u;
    ar = Avf*u;
    u = u - lambda*Ac*f + lambda^2*(ar.*(Af*f)-al.*(Ab*f)) + e*lambda*A2c*u;
    if PLOT
        plot(x,u)
        axis([min(x) max(x) -.1 1.1])
        drawnow
    end
end
end