function [c, ceq, dc, dceq] = BB_Nonlinearcon_func_LGR(x,N,D,t0,tf)

c = [];
dc = [];
dceq = [];
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);

ceq = zeros(2*N,1);
ceq(1:N,1) = D * x1' - ((tf-t0)/2)*(x2(1:N))';
ceq(N+1:2*N,1) = D * x2' -((tf-t0)/2)*(x3(1:N))';
end
