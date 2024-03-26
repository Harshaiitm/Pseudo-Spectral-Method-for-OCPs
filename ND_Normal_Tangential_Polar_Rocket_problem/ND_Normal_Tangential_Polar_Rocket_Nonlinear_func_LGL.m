function [c,ceq,dc,dceq] = ND_Normal_Tangential_Polar_Rocket_Nonlinear_func_LGL(x,M,D,problem)   %inequality constarints                      
% c = [];
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
T_max_by_W = problem.T_max_by_W;
Thrust_max = problem.Thrust_max;
q_max = problem.q_max;
a_sen_max = problem.a_sen_max;
hi = problem.hi;
hf = problem.hf;
vi = problem.vi;
vf = problem.vf;
gamma_i = problem.gamma_i;
gamma_f = problem.gamma_f;
mass_i = problem.mass_i;
theta_i = problem.theta_i;


n_length = 1/Re;
n_velocity = sqrt(Re/mu);
n_time = n_length/n_velocity;
n_mass = 1/m0;


% Decision veriables
R = x(1:M);                  % Radial position
theta = x(M+1:2*M);          % Downrange angle
V = x(2*M+1:3*M);            % tangential velocity   
gamma = x(3*M+1:4*M);        % Flight path angle
mass = x(4*M+1:5*M);         % mass
Thrust = x(5*M+1:6*M);       % Thrust
alpha = x(6*M+1:7*M);        % Angle of attack   
final_time = x(7*M+1);       % Final time

t0 = t0*n_time;
% final_time = final_time*n_time;
Isp = Isp*n_time;

h =  R/n_length - Re;
rho = rho0 * exp(-(h./h_scale));
g = mu./(R./n_length).^2;
g0 = mu/Re^2;

q = 0.5*rho.*(V/n_velocity).^2;
Drag = q.* A_ref *CD;
Drag_ref = m0*g0;
Drag = Drag./Drag_ref;


a_sen_v = (Thrust.* cos(alpha) - Drag)./(mass);
a_sen_gamma = (Thrust.* sin(alpha))./(mass);
a_sen_mag = sqrt(a_sen_v.^2 + a_sen_gamma.^2);

c = [];
% c = zeros(2*M,1);
% c(1:M,1) = (a_sen_mag.^2 - a_sen_max^2) * (mu/Re^2);
% c(M+1:2*M,1) = (q - q_max)* rho0 * (mu/Re);

ceq = zeros(5*M+8,1);
ceq(1:M,1) = D*R' - ((final_time-t0)/2)*(V.*sin(gamma))';
ceq(M+1:2*M,1) = D*theta' - ((final_time-t0)/2)*(V.*cos(gamma)./R)';
ceq(2*M+1:3*M,1) = D*V' - ((final_time-t0)/2)*((Thrust.*cos(alpha)./mass - Drag./mass) - (sin(gamma)./R.^2))';
ceq(3*M+1:4*M,1) = D*gamma' - ((final_time-t0)/2)*((Thrust.*sin(alpha)./(mass.*V)) - ((V.^2./R)-(1./R.^2)).*(cos(gamma)./V))';
ceq(4*M+1:5*M,1) = D*mass' + ((final_time-t0)/2)*(Thrust./Isp)';
ceq(5*M+1) = R(1) - (hi+Re)*n_length;
ceq(5*M+2) = R(end)-(hf+Re)*n_length;
ceq(5*M+3) = V(1)- vi*n_velocity;
ceq(5*M+4) = V(end)- vf*n_velocity;
ceq(5*M+5) = gamma(1) - gamma_i;
ceq(5*M+6) = gamma(end) - gamma_f;
ceq(5*M+7) = mass(1) - mass_i/m0;
ceq(5*M+8) = theta(1) - theta_i;

end


