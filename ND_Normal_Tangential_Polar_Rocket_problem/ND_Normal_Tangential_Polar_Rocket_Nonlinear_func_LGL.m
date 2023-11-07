function [c,ceq,dc,dceq] = ND_Normal_Tangential_Polar_Rocket_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
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
hf = problem.hf;
T_max_by_W = problem.T_max_by_W;
Thrust_max = problem.Thrust_max;
q_max = problem.q_max;
a_sen_max = problem.a_sen_max;
theta_0 = problem.theta_0;
gamma_0 = problem.gamma_0;
alpha_0 = problem.alpha_0;

n_length = 1/Re;
n_velocity = sqrt(Re/mu);
n_time = n_length/n_velocity;
n_mass = 1/m0;


% Decision veriables
R = x(1:N+1);               
theta = x(N+2:2*N+2);
V = x(2*N+3:3*N+3);
gamma = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust = x(5*N+6:6*N+6);
alpha = x(6*N+7:7*N+7);
final_time = x(7*N+8);

t0 = t0*n_time;
final_time = final_time*n_time;
Isp = Isp*n_time;

h =  R/n_length - Re;
rho = rho0 * exp(-(h./h_scale));
g = mu./(R./n_length).^2;
g0 = mu/Re;

q = 0.5*rho.*(V/n_velocity).^2;
Drag = q.* A_ref *CD;
Drag_ref = 0.5*rho_0*(mu/Re)*A_ref*CD;
Drag = Drag./Drag_ref ;

a_sen_v = (Thrust.* cos(alpha) - Drag)./(mass);
a_sen_gamma = (Thrust.* sin(alpha))./(mass);
a_sen_mag = sqrt(a_sen_v.^2 + a_sen_gamma.^2);

K = 0.5*rho0*Re^3*CD/m0;


c = zeros(2*N+2,1);
c(1:N+1,1) = (a_sen_mag.^2 - a_sen_max^2) * (mu/Re^2);
c(N+2:2*N+2,1) = (q - q_max)* 0.5 * rho0 * (mu/Re);

ceq = zeros(5*N+14,1);
ceq(1:N+1,1) = D*R' - ((final_time-t0)/2)*(V.*sin(gamma))';
ceq(N+2:2*N+2,1) = D*theta' - ((final_time-t0)/2)*(V.*cos(gamma)./R)';
ceq(2*N+3:3*N+3,1) = D*V' - ((final_time-t0)/2)*((Thrust.*cos(alpha)./mass - K*Drag./mass) - (sin(gamma)/R.^2))';
ceq(3*N+4:4*N+4,1) = D*gamma' - ((final_time-t0)/2)*((Thrust.*sin(alpha)./(mass.*V)) - ((V.^2./R)-(1/R.^2)).*(cos(gamma)./V))';
ceq(4*N+5:5*N+5,1) = D*mass' + ((final_time-t0)/2)*(Thrust./Isp)';

end


