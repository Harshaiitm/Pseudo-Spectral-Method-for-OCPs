function Objective = Objective_LGR(x,M,weights,t0,tf,t)
x1 = x(2:M);
x2 = x(M+2:2*M);
x3 = x(2*M+2:3*M);
xr = 5*sin(t(2:M))';
vr = 5*cos(t(2:M))';

Objective = ((tf-t0)/2)*(sum(((x1-xr).*(x1-xr))'.*weights + ((x2-vr).*(x2-vr))'.*weights + (x3.*x3)'.*weights*0.0001));
end