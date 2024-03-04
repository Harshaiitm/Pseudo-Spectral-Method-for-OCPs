function [c, ceq, dc, dceq] = Nonlinearcon_LGR(x,x0,M,D,t0,tf)
c = [];
dc = [];
dceq = [];
x1 = x(1:M);
x2 = x(M+1:2*M);
x3 = x(2*M+1:3*M);

ceq = zeros(2*M,1);
ceq(1:M,1) = D*[x0(1) x1]' - ((tf-t0)/2) * x2';
ceq(M+1:2*M,1) = D*[x0(M+1) x2]' -((tf-t0)/2) * x3';
end
