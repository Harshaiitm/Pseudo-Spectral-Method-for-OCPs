function [c,ceq,dc,dceq] = ND_Radial_Transverse_Polar_Rocket_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
dc = [];
dceq = [];

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

n_length = 1/Re;
n_velocity = sqrt(Re/mu);
n_time = n_length/n_velocity;
n_mass = 1/m0;


% Decision veriables
R = x(1:N+1);
theta = x(N+2:2*N+2);
V_r = x(2*N+3:3*N+3);
V_theta = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust_r = x(5*N+6:6*N+6);
Thrust_theta = x(6*N+7:7*N+7);
final_time = x(7*N+8);

t0 = t0*n_time;
% final_time = final_time*n_time;
Isp = Isp*n_time;

h =  R/n_length - Re;
rho = rho0 * exp(-(h./h_scale));
g = mu./(R./n_length).^2;
g0 = mu/Re^2;
% q = 0.5*rho.*((V_r/n_velocity).^2 + (V_theta./n_velocity).^2);
Thrust_mag = (Thrust_r.^2 + Thrust_theta.^2).^0.5;
Drag_r = - 0.5*rho.* V_r.*(V_r.^2 + V_theta.^2).^0.5 *A_ref *CD;
Drag_ref = m0*g0;
Drag_r = Drag_r./Drag_ref;
Drag_theta = - 0.5*rho.* V_theta.*(V_r.^2 + V_theta.^2).^0.5 *A_ref *CD;
Drag_theta = Drag_theta./Drag_ref;

a_sen_r = (Thrust_r + Drag_r)./(mass);
a_sen_theta = (Thrust_theta + Drag_theta)./(mass);
a_sen_mag = sqrt(a_sen_r.^2 + a_sen_theta.^2);

c = [];
% c = zeros(2*N+2,1);
% c(1:N+1,1) = (a_sen_mag.^2 - a_sen_max^2) * (mu/Re^2);
% c(N+2:2*N+2,1) = (q - q_max)* 0.5 * rho0 * (mu/Re);

ceq = zeros(5*N+8,1);
ceq(1:N+1,1) = D*R' - ((final_time-t0)/2)*(V_r)';
ceq(N+2:2*N+2,1) = D*theta' - ((final_time-t0)/2)*(V_theta./R)';
ceq(2*N+3:3*N+3,1) = D*V_r' - ((final_time-t0)/2)*((Thrust_r./mass + Drag_r./mass) - (1./R.^2)+ ((V_theta.^2)./R))';
ceq(3*N+4:4*N+4,1) = D*V_theta' - ((final_time-t0)/2)*((Thrust_r./mass + Drag_r./mass) - (1./R.^2)+ (V_r.*V_theta)./R)';
ceq(4*N+5:5*N+5,1) = D*mass' + ((final_time-t0)/2)*(((Thrust_r.^2+Thrust_theta.^2).^(0.5))./Isp)';
% ceq(5*N+6) = R(end)-(hf+Re)*n_length;
% ceq(5*N+7) = V(end)- sqrt(mu/(hf+Re))*n_velocity;
% ceq(5*N+8) = gamma(end)-0;


end


