function [x,u] = Godunov(H,K,POS,PLOT)
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
u2 = u;

for i = 1:(1/k)
    for j = 2:xl-1
        u2(j) = u(j) - 1/2*lambda*((Reimann(0,u(j),u(j+1)))^2-(Reimann(0,u(j-1),u(j)))^2);
    end
    u = u2;
    if PLOT
        plot(x,u)
        axis([min(x) max(x) 0 1])
        drawnow
    end
end
end