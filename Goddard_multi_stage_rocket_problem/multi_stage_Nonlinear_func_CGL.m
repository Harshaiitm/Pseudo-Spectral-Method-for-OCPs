function [c,ceq,dc,dceq] = multi_stage_Nonlinear_func_CGL(x,M,D,problem)   %inequality constarints                      

c = [];
dc = [];
dceq = [];

% Decision veriables
h_1 = x(1:M);
v_1 = x(M+1:2*M);
mass_1 = x(2*M+1:3*M);
Thrust_1 = x(3*M+1:4*M);
h_2 = x(4*M+1:5*M);
v_2 = x(5*M+1:6*M);
mass_2 = x(6*M+1:7*M);
Thrust_2 = x(7*M+1:8*M);
stage_time = x(8*M+1);
final_time = x(8*M+2);

Re = problem.Re;
h_scale = problem.h_scale;
mu =problem.mu;
m0 = problem.m0;
m0_1 = problem.m0_1;
m0_2 = problem.m0_2;
A_ref = problem.A_ref;
CD = problem.CD;
rho0 = problem.rho0;
g0 = problem.g0;
Isp = problem.Isp;
t0 = problem.t0;
ts = problem.ts;
tf = problem.tf;

r_1 = h_1 + Re;
r_2 = h_2 + Re;
rho_1 = rho0 * exp(-(1/h_scale).*(h_1));
rho_2 = rho0 * exp(-(1/h_scale).*(h_2));
g_1 = mu./r_1.^2;
g_2 = mu./r_2.^2;
Drag_1 = 0.5*rho_1.* v_1.^2 *A_ref *CD;
Drag_2 = 0.5*rho_2.* v_2.^2 *A_ref *CD;


ceq = zeros(6*M+8,1);
ceq(1:M,1) = D*h_1'-((stage_time-t0)/2)* (v_1)';
ceq(M+1:2*M,1)= D*v_1' - ((stage_time-t0)/2)*((Thrust_1-Drag_1)./mass_1 - mu./r_1.^2)';
ceq(2*M+1:3*M,1) = D*mass_1'+((stage_time-t0)/2)*(Thrust_1./(g0.*Isp))';
ceq(3*M+1:4*M,1) = D*h_2'-((final_time-stage_time)/2)* (v_2)';
ceq(4*M+1:5*M,1)= D*v_2' - ((final_time-stage_time)/2)*((Thrust_2-Drag_2)./mass_2 - mu./r_2.^2)';
ceq(5*M+1:6*M,1) = D*mass_2'+((final_time-stage_time)/2)*(Thrust_2./(g0.*Isp))';
ceq(6*M+1) = 0-h_1(1);
ceq(6*M+2) = h_2(1)-h_1(end);
ceq(6*M+3) = 0-v_1(1);
ceq(6*M+4) = v_2(1)-v_1(end);
ceq(6*M+5) = mass_1(1)-5000;
ceq(6*M+6) = mass_1(end)-3200;
ceq(6*M+7) = mass_2(1)-2000;
ceq(6*M+8) = mass_2(end)-800;
 
end




