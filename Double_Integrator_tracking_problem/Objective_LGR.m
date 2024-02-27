function Objective = Objective_LGR(x,x0,N,weights,t0,tf,t)
x1 = [x0(1) x(2:N+1)];
x2 = [x0(N+2) x(N+3:2*N+2)];
x3 = [NaN x(2*N+4:3*N+3)];
xr = 5*sin(t');
v = 5*cos(t');

Objective = ((tf-t0)/2)*(sum(((x1(2:N+1)-xr(2:N+1)).*(x1(2:N+1)-xr(2:N+1)))'.*weights + ((x2(2:N+1)-v(2:N+1)).*(x2(2:N+1)-v(2:N+1)))'.*weights + (x3(2:N+1).*x3(2:N+1))'.*weights*0.0001));
end