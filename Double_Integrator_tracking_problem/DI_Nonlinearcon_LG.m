function [c, ceq, dc, dceq] = DI_Nonlinearcon_LG(x,M,D,problem)
c = [];
dc = [];
dceq = [];

xi = problem.xi;
vi = problem.vi;
t0 = problem.t0;
tf = problem.tf;

x1 = x(1:M);
x2 = x(M+1:2*M);
x3 = x(2*M+1:3*M);

ceq = zeros(2*M,1);
ceq(1:M,1) = D*[xi x1]' - ((tf-t0)/2) * x2';
ceq(M+1:2*M,1) = D*[vi x2]' -((tf-t0)/2) * x3';

end
