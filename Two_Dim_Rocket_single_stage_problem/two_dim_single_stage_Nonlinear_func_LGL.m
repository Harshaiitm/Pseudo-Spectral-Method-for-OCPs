function [c,ceq,dc,dceq] = two_dim_single_stage_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
% c = [];
dc = [];
dceq = [];

% Decision veriables
R_x = x(1:N+1);
R_y = x(N+2:2*N+2);
V_x = x(2*N+3:3*N+3);
V_y = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust_x = x(5*N+6:6*N+6);
Thrust_y = x(6*N+7:7*N+7);
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


R_mag = (R_x.^2+R_y.^2).^0.5;
h =  R_mag - Re;
V_mag = (V_x.^2+V_y.^2).^0.5;
rho = rho0 * exp(-(1/h_scale).*(h));
Thrust_mag = (Thrust_x.^2+Thrust_y.^2).^0.5;
g_x = -mu*(R_x)./((R_x).^2+(R_y).^2).^(3/2);
g_y = -mu*(R_y)./((R_x).^2+(R_y).^2).^(3/2);
Drag_x = - 0.5*rho.* V_x.*(V_x.^2+V_y.^2).^0.5 *A_ref *CD;
Drag_y = - 0.5*rho.* V_y.*(V_x.^2+V_y.^2).^0.5 *A_ref *CD;
q = 0.5*rho.*(V_x.^2+V_y.^2);
a_sen_x = (Thrust_x+Drag_x)./mass;
a_sen_y = (Thrust_y+Drag_y)./mass;
a_sen_mag = (a_sen_x.^2+a_sen_y.^2).^0.5;


c = zeros(3*N+3,1);
c(1:N+1,1) = Thrust_mag.^2 - Thrust_max^2;
c(N+2:2*N+2,1) = a_sen_mag.^2 - a_sen_max^2;
c(2*N+3:3*N+3,1) = q - q_max;

ceq = zeros(5*N+15,1);
ceq(1:N+1,1) = D*R_x' - ((final_time-t0)/2)*(V_x)';
ceq(N+2:2*N+2,1) = D*R_y' - ((final_time-t0)/2)*(V_y)';
ceq(2*N+3:3*N+3,1) = D*V_x' - ((final_time-t0)/2)*((Thrust_x+Drag_x)./mass + g_x)';
ceq(3*N+4:4*N+4,1) = D*V_y' - ((final_time-t0)/2)*((Thrust_y+Drag_y)./mass + g_y)';
ceq(4*N+5:5*N+5,1) = D*mass' + ((final_time-t0)/2)*((Thrust_x.^2+Thrust_y.^2).^0.5./(g0.*Isp))';
ceq(5*N+6,1) = Re*sqrt(1/2) - R_x(1);
ceq(5*N+7,1) = Re*sqrt(1/2) -  R_y(1);
ceq(5*N+8,1) = R_mag(end)^2-(R_x(end)^2+R_y(end)^2);
ceq(5*N+9,1) = (hf+Re)^2-(R_x(end)^2+R_y(end)^2);
ceq(5*N+10,1) = 0-V_x(1);
ceq(5*N+11,1) = 0-V_y(1);
ceq(5*N+12,1) = V_mag(end)^0.5 -(V_x(end)^2+V_y(end)^2)^0.5;
ceq(5*N+13,1) = (mu/(hf+Re))^0.5 -(V_x(end)^2+V_y(end)^2)^0.5;
ceq(5*N+14,1) = 0-(V_x(end)*R_x(end)+V_y(end)*R_y(end));
ceq(5*N+15,1) = m0-mass(1);
end



