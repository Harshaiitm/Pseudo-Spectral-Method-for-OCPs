function [c,ceq,dc,dceq] = single_stage_Nonlinear_func_LGR(x,M,D,problem)                         
c = [];
dc = [];
dceq = [];

% Decision veriables
h = x(1:M);
v = x(M+1:2*M);
mass = x(2*M+1:3*M);
Thrust = x(3*M+1:4*M);
final_time = x(4*M+1);

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
hi = problem.hi;
vi = problem.vi;
mass_i = problem.mass_i;
mass_f = problem.mass_f;
Thrust_i = problem.Thrust_i;
x0 = problem.x0;

r = h + Re;
rho = rho0 * exp(-(1/h_scale).*(h));
g = mu./r.^2;
Drag = 0.5*rho.* v.^2 *A_ref *CD;
% mp = m0-mass;


ceq = zeros(3*M+1,1);
ceq(1:M,1) = D*[hi h]'-((final_time-t0)/2)* v';
ceq(M+1:2*M,1)=D*[vi v]' - ((final_time-t0)/2)*((Thrust - Drag)./mass - mu./r.^2)';
ceq(2*M+1:3*M,1) = D*[mass_i mass]'+((final_time-t0)/2)*(Thrust./(g0.*Isp))';
ceq(3*M+1) = mass(end)-mass_f;
end




