function Objective = Objective_LG(x,N,weights,t0,tf,t)
x1 = x(1:N);                   % position  
x1(N+1)=x1(1)+weights(1:N)'*(x1(1:N))';
x2 = x(N+2:2*N+1);               % velocity 
x2(N+1)=x2(1)+weights(1:N)'*(x2(1:N))';
x3 = x(2*N+3:3*N+3);             % accleration
x3(N+1)=x3(1)+weights(1:N)'*(x3(1:N))';
xr = 5*sin(t');
v = 5*cos(t');

Objective = ((tf-t0)/2)*(sum(((x1(1:N)-xr(1:N)).*(x1(1:N)-xr(1:N)))'.*weights+((x2(1:N)-v(1:N)).*(x2(1:N)-v(1:N)))'.*weights+(x3(1:N).*x3(1:N))'.*weights*0.0001));
end