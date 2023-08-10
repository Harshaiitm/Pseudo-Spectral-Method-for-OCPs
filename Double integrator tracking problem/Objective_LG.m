function Objective = Objective_LG(x,N,weights,t)
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
xr = 5*sin(t');
v = 5*cos(t');

Objective = sum(((x1(1:N)-xr(1:N)).*(x1(1:N)-xr(1:N)))'.*weights+((x2(1:N)-v(1:N)).*(x2(1:N)-v(1:N)))'.*weights+(x3(1:N).*x3(1:N))'.*weights*0.0001);
end