function [c,ceq,dc,dceq] = single_stage_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
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

ceq = zeros(3*N+8,1);
ceq(1:N+1,1) = D*h-((final_time-t0)/2)* (v);
ceq(N+2:2*N+2,1)=D*v - ((final_time-t0)/2)*((Thrust-Drag)./mass - mu./r.^2);
ceq(2*N+3:3*N+3,1) = D*mass+((final_time-t0)/2)*(Thrust./(g0.*Isp));
ceq(3*N+4) = 0-h(1);
ceq(3*N+5) =0-v(1);
ceq(3*N+6) =m0-mass(1);
ceq(3*N+7) = mass(end)-2000;
ceq(3*N+8) = Thrust(1)-m0*g0*2;

end



