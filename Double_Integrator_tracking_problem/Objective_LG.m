function Objective = Objective_LG(x,M,weights,t0,tf,t)
x1 = x(1:M);                   % position  
x2 = x(M+1:2*M);               % velocity 
x3 = x(2*M+1:3*M);             % accleration
xr = 5*sin(t');
vr = 5*cos(t');

Objective = ((tf-t0)/2)*(sum(((x1-xr).*(x1-xr))'.*weights + ((x2-vr).*(x2-vr))'.*weights + (x3.*x3)'.*weights*0.0001));
end