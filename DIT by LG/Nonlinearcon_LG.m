function [c, ceq, dc, dceq] = Nonlinearcon_LG(x,N,D,t0,tf,weights)
c = [];
dc = [];
dceq = [];
x1 = x(1:N+1);
% x1(N+2) = x1(1)+((tf-t0)/2)*(x1(1:N))*weights;
x2 = x(N+2:2*N+2);
% x2(N+2) = x2(1)+((tf-t0)/2)*(x2(1:N))*weights;
x3 = x(2*N+3:3*N+3);
% x3(N+2) = x3(1)+((tf-t0)/2)*(x3(1:N))*weights;
ceq = zeros(2*N,1);
ceq(1:N,1) = D * x1' - ((tf-t0)/2)*(x2(1:N))';
ceq(N+1:2*N,1) = D * x2' -((tf-t0)/2)*(x3(1:N))'; 
ceq(2*N+1) = x(1)-0;            % starting point
ceq(2*N+2) = x(N+2)-5;
end