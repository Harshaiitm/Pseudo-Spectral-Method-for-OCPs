function [c,ceq,dc,dceq] = single_stage_Nonlinear_func_LGR(x,N,D,problem)   %inequality constarints                      
c = [];
dc = [];
dceq = [];

% Decision veriables
h = x(1:N+1);
v = x(N+2:2*N+2);
mass = x(2*N+3:3*N+3);
Thrust = x(3*N+4:4*N+4);
final_time = x(4*N+5);

Re = problem.Re;
h_scale = problem.h_scale;
mu =problem.mu;
m0 = problem.m0;
A_ref = problem.A_ref;
CD = problem.CD;
rho0 = problem.rho0;
g0 = problem.g0;
Isp = problem.Isp;
t0 = problem.t0;
tf = problem.tf;

r = h + Re;
rho = rho0 * exp(-(1/h_scale).*(h));
g = mu./r.^2;
Drag = 0.5*rho.* v.^2 *A_ref *CD;
% mp = m0-mass;

ceq = zeros(3*N+5,1);
ceq(1:N,1) = D*h'-((final_time-t0)/2)* (v(1:N))';
ceq(N+1:2*N,1)=D*v' - ((final_time-t0)/2)*((Thrust(1:N)-Drag(1:N))./mass(1:N) - mu./r(1:N).^2)';
ceq(2*N+1:3*N,1) = D*mass'+((final_time-t0)/2)*(Thrust(1:N)./(g0.*Isp))';
ceq(3*N+1) = 0-h(1);
ceq(3*N+2) =0-v(1);
ceq(3*N+3) =m0-mass(1);
ceq(3*N+4) = mass(end)-2000;
ceq(3*N+5) = Thrust(1)-m0*g0*2;

end




