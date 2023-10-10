function [c,ceq,dc,dceq] = two_dim_single_stage_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
% c = [];
dc = [];
dceq = [];

% Decision veriables
h_x = x(1:N+1);
h_y = x(N+2:2*N+2);
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
tf = problem.tf;
hf = problem.hf;
T_max_by_W = problem.T_max_by_W;
Thrust_max = problem.Thrust_max;
q_max = problem.q_max;
a_sen_max = problem.a_sen_max;


h = (h_x.^2+h_y.^2).^0.5;
R = Re + h;
V = (V_x.^2+V_y.^2).^0.5;
rho = rho0 * exp(-(1/h_scale).*(R));
Thrust = (Thrust_x.^2+Thrust_y.^2).^0.5;
g_x = -mu*(h_x+Re)./((h_x+Re).^2+(h_y+Re).^2).^(3/2);
g_y = -mu*(h_y+Re)./((h_x+Re).^2+(h_y+Re).^2).^(3/2);
Drag_x = - 0.5*rho.* V_x.*(V_x.^2+V_y.^2).^2 *A_ref *CD;
Drag_y = - 0.5*rho.* V_y.*(V_x.^2+V_y.^2).^2 *A_ref *CD;
q = 0.5*rho.*(V_x.^2+V_y.^2);
a_sen_x = (Thrust_x+Drag_x)./mass;
a_sen_y = (Thrust_y+Drag_y)./mass;
a_sen = (a_sen_x.^2+a_sen_y.^2).^0.5;


c = zeros(3*N+3,1);
c(1:N+1,1) = Thrust.^2 - Thrust_max^2;
c(N+2:2*N+2,1) = a_sen.^2 - a_sen_max^2;
c(2*N+3:3*N+3,1) = q - q_max;

ceq = zeros(6*N+14,1,1);
ceq(1:N+1,1) = D*h_x' - ((final_time-t0)/2)*(V_x)';
ceq(N+2:2*N+2,1) = D*h_y' - ((final_time-t0)/2)*(V_y)';
ceq(2*N+3:3*N+3,1)=D*V_x' - ((final_time-t0)/2)*((Thrust_x-Drag_x)./mass - g_x)';
ceq(3*N+4:4*N+4,1)=D*V_y' - ((final_time-t0)/2)*((Thrust_y-Drag_y)./mass - g_y)';
ceq(4*N+5:5*N+5,1) = D*mass' + ((final_time-t0)/2)*((Thrust_x.^2+Thrust_y.^2).^0.5./(g0.*Isp))';
ceq(6*N+7,1) = 0-h_x(1);
ceq(6*N+8,1) = 0-h_y(1);
ceq(6*N+9,1) = 0-h(1);
ceq(6*N+10,1) = hf - h(end);
% ceq(6*N+10,1) = q - q_max;
% ceq(6*N+11,1) = a_sen^2 - a_sen_max^2;
ceq(6*N+11,1) = h(end)-(h_x(end)^2+h_y(end)^2);
ceq(6*N+12,1) = 0-V_x(1);
ceq(6*N+13,1) = 0-V_y(1);
ceq(6*N+14,1) = 0-V(1);
ceq(6*N+15,1) = V(end)-(V_x(end)^2+V_y(end)^2);
ceq(6*N+16,1) = 0-(V_x(end)*h_x(end)+V_y(end)*h_y(end));
ceq(6*N+17,1) = m0-mass(1);
ceq(6*N+18,1) = mass(end)-1000;
ceq(6*N+19,1) = (mu/R(end))^0.5 - V(end);
% ceq(6*N+19,1) = (Thrust_x.^2+Thrust_y.^2).^0.5-Thrust;
% ceq(6*N+20,1) = Thrust_x(1)-m0*g0*T_max_by_W; 
% ceq(6*N+21,1) = Thrust_y(1)-m0*g0*T_max_by_W;
end




