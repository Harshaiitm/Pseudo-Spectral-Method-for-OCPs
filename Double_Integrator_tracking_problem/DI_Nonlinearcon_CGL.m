function [c, ceq, dc, dceq] = DI_Nonlinearcon_CGL(x,M,D,problem)
c = [];                         %inequality constarints
dc = [];
dceq = [];

xi = problem.xi;
vi = problem.vi;
t0 = problem.t0;
tf = problem.tf;

x1 = x(1:M);
x2 = x(M+1:2*M);
x3 = x(2*M+1:3*M);

ceq = zeros(2*M+2,1);              % equality constarints
ceq(1:M,1) = D * x1' - ((tf-t0)/2)*x2';
ceq(M+1:2*M,1) = D * x2' -((tf-t0)/2)*x3';
ceq(2*M+1) = xi - x1(1);                 % initial conditions
ceq(2*M+2) = vi - x2(1);
end