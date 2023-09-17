function [c, ceq, dc, dceq] = Nonlinearcon_LGR(x,N,D,t0,tf)

c = [];
dc = [];
dceq = [];
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);

c = zeros(1,N+1);
c(1,1:N+1) = 0-x1;

ceq = zeros(1,2*N+6);
ceq(1,1:N) = D * x1' - ((tf-t0)/2)*x2(1:N)';
ceq(1,N+1:2*N) = D * x2' -((tf-t0)/2)*x3(1:N)';
ceq(2*N+3) = x1(end)-300;
ceq(2*N+4) = x2(1)-0;
ceq(2*N+5) = x2(end)-0;
ceq(2*N+6) =x1(end)-(x3(end)/2)*(x4^2);
end
