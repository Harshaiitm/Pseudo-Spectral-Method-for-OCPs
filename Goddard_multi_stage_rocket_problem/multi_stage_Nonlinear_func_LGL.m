function [c,ceq,dc,dceq] = multi_stage_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
dc = [];
dceq = [];

% Decision veriables
h_1 = x(1:N+1);
v_1 = x(N+2:2*N+2);
mass_1 = x(2*N+3:3*N+3);
Thrust_1 = x(3*N+4:4*N+4);
h_2 = x(4*N+5:5*N+5);
v_2 = x(5*N+6:6*N+6);
mass_2 = x(6*N+7:7*N+7);
Thrust_2 = x(7*N+8:8*N+8);
stage_time = x(8*N+9);
final_time = x(8*N+10);

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


ceq = zeros(6*N+16,1);
ceq(1:N+1,1) = D*h_1'-((stage_time-t0)/2)* (v_1)';
ceq(N+2:2*N+2,1)=D*v_1' - ((stage_time-t0)/2)*((Thrust_1-Drag_1)./mass_1 - mu./r_1.^2)';
ceq(2*N+3:3*N+3,1) = D*mass_1'+((stage_time-t0)/2)*(Thrust_1./(g0.*Isp))';
ceq(3*N+4:4*N+4,1) = D*h_2'-((final_time-stage_time)/2)* (v_2)';
ceq(4*N+5:5*N+5,1)=D*v_2' - ((final_time-stage_time)/2)*((Thrust_2-Drag_2)./mass_2 - mu./r_2.^2)';
ceq(5*N+6:6*N+6,1) = D*mass_2'+((final_time-stage_time)/2)*(Thrust_2./(g0.*Isp))';
ceq(6*N+7) = 0-h_1(1);
ceq(6*N+8) = h_2(1)-h_1(end);
ceq(6*N+9) = 0-v_1(1);
ceq(6*N+10) = v_2(1)-v_1(end);
ceq(6*N+11) = mass_1(1)-5000;
ceq(6*N+12) = mass_1(end)-3800;
ceq(6*N+13) = Thrust_1(1)-m0*g0*2;
ceq(6*N+14) = mass_2(1)-2000;
ceq(6*N+15) = mass_2(end)-800;
ceq(6*N+16) = Thrust_2(1)-m0_2*g0*2;
 
% c = [];

c = zeros(2,1);
c(1) =mass_2(1)-mass_1(end);
c(2) = stage_time-final_time+0.100;
 
end




