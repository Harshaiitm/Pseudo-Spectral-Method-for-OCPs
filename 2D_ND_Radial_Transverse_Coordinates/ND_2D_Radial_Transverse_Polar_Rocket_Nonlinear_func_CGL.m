function [c,ceq,dc,dceq] = ND_2D_Radial_Transverse_Polar_Rocket_Nonlinear_func_CGL(x,M,D,problem)   %inequality constarints                      
dc = [];
dceq = [];


Omega_z = 2*pi/(24*60*60);

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
Vi = problem.Vi;
theta_i = problem.theta_i;
% Omega_z = problem.Omega_z;
hf = problem.hf;
Vf = problem.Vf;
T_max_by_W = problem.T_max_by_W;
Thrust_max = problem.Thrust_max;
q_max = problem.q_max;
a_sen_max = problem.a_sen_max;

n_length = 1/Re;
n_velocity = sqrt(Re/mu);
n_time = n_length/n_velocity;
n_mass = 1/m0;
n_thrust = 1/(m0*g0);

% Decision veriables
R = x(0*M+1:1*M);               
theta = x(1*M+1:2*M);
V_r = x(2*M+1:3*M);
V_theta = x(3*M+1:4*M);
mass = x(4*M+1:5*M);
Thrust_r = x(5*M+1:6*M);
Thrust_theta = x(6*M+1:7*M);
final_time = x(7*M+1);

t0 = t0*n_time;
% final_time = final_time*n_time;
Isp = Isp*n_time;

h =  R/n_length - Re;
rho = rho0 * exp(-(h./h_scale));
V_mag = sqrt((V_r/n_velocity).^2 + (V_theta/n_velocity).^2);
q = 0.5*rho.* V_mag.^2;
Thrust_mag = sqrt((Thrust_r/n_thrust).^2 + (Thrust_theta/n_thrust).^2);
Drag_r = - 0.5*rho.*(V_r/n_velocity).*V_mag*A_ref *CD;
n_drag = 1/(m0*g0);
Drag_r = Drag_r.* n_drag;
Drag_theta = - 0.5*rho.* (V_theta/n_velocity).*V_mag*A_ref *CD;
Drag_theta = Drag_theta.* n_drag;

a_sen_r = (Thrust_r/n_thrust + Drag_r/n_drag)./(mass/n_mass);
a_sen_theta = (Thrust_theta/n_thrust + Drag_theta/n_drag)./(mass/n_mass);
a_sen_mag = sqrt(a_sen_r.^2 + a_sen_theta.^2);

% c = [];
c = zeros(3*M,1);
c(1:M,1) = (a_sen_mag.^2 - a_sen_max^2);
c(M+1:2*M,1) = (q - q_max);
c(2*M+1:3*M,1) = (Thrust_r.^2 + Thrust_theta.^2) -Thrust_max^2;

ceq = zeros(5*M+11,1);
ceq(1:M,1) = D*R' - ((final_time-t0)/2)*(V_r)';
ceq(M+1:2*M,1) = D*theta' - ((final_time-t0)/2)*(V_theta./R)';
ceq(2*M+1:3*M,1) = D*V_r' - ((final_time-t0)/2)*((Thrust_r./mass + Drag_r./mass) - (1./R.^2)+ ((V_theta.^2)./R))';
ceq(3*M+1:4*M,1) = D*V_theta' - ((final_time-t0)/2)*((Thrust_theta./mass + Drag_theta./mass) - (V_r.*V_theta)./R)';
ceq(4*M+1:5*M,1) = D*mass' + ((final_time-t0)/2)*(sqrt(Thrust_r.^2 + Thrust_theta.^2)/Isp)';

ceq(5*M+2) = R(1) - (hi+Re)*n_length;
ceq(5*M+3) = R(end) - (hf+Re)*n_length;
ceq(5*M+4) = V_r(1) - (Vi)*n_velocity;
ceq(5*M+5) = V_theta(1) - 0;
ceq(5*M+6) = V_r(end) - 0;
ceq(5*M+7) = V_theta(end) - Vf * n_velocity;
ceq(5*M+8) = mass(1) - m0*n_mass;
ceq(5*M+9) = theta(1) - theta_i;
% 
ceq(5*M+10) = D(1,:)*theta'- Omega_z; 
% ceq(5*M+11) = D(M,:)*theta' - Vf*n_velocity/(Re*n_length + hf*n_length);

end


