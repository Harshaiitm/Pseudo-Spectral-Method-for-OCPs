function Objective = Objective_CGL(x,N,weights,t)
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
xr = 5*sin(t');  % referece for position to follow 
v = 5*cos(t');   % referece for velocity to follow   

Objective = sum(((x1-xr).*(x1-xr))'.*weights+((x2-v).*(x2-v))'.*weights+(x3.*x3)'.*weights*0.0001);
end 