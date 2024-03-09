function [c,ceq,dc,dceq] = single_stage_Nonlinear_func_LGR(x,x0,M,D,problem)   %inequality constarints                      
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
tf = problem.tf;

r = h + Re;
rho = rho0 * exp(-(1/h_scale).*(h));
g = mu./r.^2;
Drag = 0.5*rho.* v.^2 *A_ref *CD;
% mp = m0-mass;


ceq = zeros(3*M+5,1);
ceq(1:M,1) = D*[x0(1) h]'-((final_time-t0)/2)* v';
ceq(M+1:2*M,1)=D*[x0(M+1) v]' - ((final_time-t0)/2)*((Thrust - Drag)./mass - mu./r.^2)';
ceq(2*M+1:3*M,1) = D*[x0(2*M+1) mass]'+((final_time-t0)/2)*(Thrust./(g0.*Isp))';
ceq(3*M+1) = 0-h(1);
ceq(3*M+2) =0-v(1);
ceq(3*M+3) =m0-mass(1);
ceq(3*M+4) = mass(end)-2000;
ceq(3*M+5) = Thrust(1)-m0*g0*2;
end




