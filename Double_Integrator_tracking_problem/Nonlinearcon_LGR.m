function [c, ceq, dc, dceq] = Nonlinearcon_LGR(x,x0,N,D,t0,tf)
c = [];
dc = [];
dceq = [];
x1 = [x0(1) x(2:N+1)];
x2 = [x0(N+2) x(N+3:2*N+2)];
x3 = [x0(2*N+3) x(2*N+4:3*N+3)];
ceq = zeros(2*N+2,1);
ceq(1:N,1) = D * x1' - ((tf-t0)/2)*(x2(2:N+1))';
ceq(N+1:2*N,1) = D * x2' -((tf-t0)/2)*(x3(2:N+1))';
ceq(2*N+1) = x(1)-0;            % starting point
ceq(2*N+2) = x(N+2)-5;
end
