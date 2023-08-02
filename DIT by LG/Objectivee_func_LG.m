function Objective = Objectivee_func_LG(x,N,weights,t,t0,tf)
x1 = x(1:N+1);
% x1(N+2) = x1(1)+((tf-t0)/2)*(x1(1:N))*weights;
x2 = x(N+2:2*N+2);
% x2(N+2) = x2(1)+((tf-t0)/2)*(x2(1:N))*weights;
x3 = x(2*N+3:3*N+3);
% x3(N+2) = x3(1)+((tf-t0)/2)*(x3(1:N))*weights;
xr = 5*sin(t');
v = 5*cos(t');

Objective = @(x) sum(((x1(1:N)-xr(1:N)).*(x1(1:N)-xr(1:N)))'.*weights+((x2(1:N)-v(1:N)).*(x2(1:N)-v(1:N)))'.*weights+(x3(1:N).*x3(1:N))'.*weights);
end