function [c,ceq,dc,dceq] = Radial_Transverse_Polar_Rocket_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
% c = [];
dc = [];
dceq = [];

% Decision veriables
R = x(1:N+1);
theta = x(N+2:2*N+2);
V_r = x(2*N+3:3*N+3);
V_theta = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust_r = x(5*N+6:6*N+6);
Thrust_theta = x(6*N+7:7*N+7);
final_time = x(7*N+8);

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
hf = problem.hf;
T_max_by_W = problem.T_max_by_W;
Thrust_max = problem.Thrust_max;
q_max = problem.q_max;
a_sen_max = problem.a_sen_max;


h =  R - Re;
rho = rho0 * exp(-(1/h_scale).*(h));
Thrust_mag = (Thrust_r.^2 + Thrust_theta.^2).^0.5;
g = mu./(R).^2;
Drag_r = - 0.5*rho.* V_r.*(V_r.^2 + V_theta.^2).^0.5 *A_ref *CD;
Drag_theta = - 0.5*rho.* V_theta.*(V_r.^2 + V_theta.^2).^0.5 *A_ref *CD;
q = 0.5*rho.*(V_r.^2 + V_theta.^2);
a_sen_r = (Thrust_r + Drag_r)./mass;
a_sen_theta = (Thrust_theta + Drag_theta)./mass;
a_sen_mag = (a_sen_r.^2 + a_sen_theta.^2).^0.5;


c = zeros(3*N+3,1);
c(1:N+1,1) = Thrust_mag.^2 - Thrust_max^2;
c(N+2:2*N+2,1) = a_sen_mag.^2 - a_sen_max^2;
c(2*N+3:3*N+3,1) = q - q_max;

ceq = zeros(5*N+11,1);
ceq(1:N+1,1) = D*R' - ((final_time-t0)/2)*(V_r)';
ceq(N+2:2*N+2,1) = D*theta' - ((final_time-t0)/2)*(V_theta/R)';
ceq(2*N+3:3*N+3,1) = D*V_r' - ((final_time-t0)/2)*((Thrust_r./mass + Drag_r./mass) + g - V_theta.^2./R)';
ceq(3*N+4:4*N+4,1) = D*V_theta' - ((final_time-t0)/2)*((Thrust_theta./mass + Drag_theta./mass) + (V_r.*V_theta)./R)';
ceq(4*N+5:5*N+5,1) = D*mass' + ((final_time-t0)/2)*((Thrust_r.^2+Thrust_theta.^2).^0.5./(g0.*Isp))';
ceq(5*N+6,1) = Re - R(1);
ceq(5*N+7,1) = (hf+Re) -  R(end);
ceq(5*N+8,1) = 0 - V_r(1);
ceq(5*N+9,1) = 0 - V_r(end);
ceq(5*N+10,1) = (mu/(hf+Re))^0.5- V_theta(end);
ceq(5*N+11,1) = m0-mass(1);
end


