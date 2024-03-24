function [c,ceq,dc,dceq] = Normal_Tangential_Polar_Rocket_Nonlinear_func_LGL(x,M,D,problem)   %inequality constarints                      


% c = [];
dc = [];
dceq = [];

% Decision veriables
R = x(1:M);               
theta = x(M+1:2*M);
V = x(2*M+1:3*M);
gamma = x(3*M+1:4*M);
mass = x(4*M+1:5*M);
Thrust = x(5*M+1:6*M);
alpha = x(6*M+1:7*M);
final_time = x(7*M+1);


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
Ri = problem.Ri;
Vi = problem.Vi;
Vf = problem.Vf;
theta_i = problem.theta_i;
gamma_i = problem.gamma_i;
gamma_f = problem.gamma_f;
alpha_i = problem.alpha_i;
mass_i = problem.mass_i; 
T_max_by_W = problem.T_max_by_W;
Thrust_max = problem.Thrust_max;
q_max = problem.q_max;
a_sen_max = problem.a_sen_max;


h =  R - Re;
rho = rho0 * exp(-(1/h_scale).*(h));
g = mu./(R).^2;
q = 0.5*rho.*(V.^2);
Drag = q.* A_ref *CD;
a_sen_v = (Thrust.* cos(alpha) - Drag)./mass;
a_sen_gamma = (Thrust.* sin(alpha))./mass;
a_sen_mag = sqrt(a_sen_v.^2 + a_sen_gamma.^2);

c = zeros(2*M,1);
c(1:M,1) = a_sen_mag.^2 - a_sen_max^2;
c(M+1:2*M,1) = q - q_max;

ceq = zeros(5*M+9,1);
ceq(1:M,1) = D*R' - ((final_time-t0)/2)*(V.*sin(gamma))';
ceq(M+1:2*M,1) = D*theta' - ((final_time-t0)/2)*(V.*cos(gamma)./R)';
ceq(2*M+1:3*M,1) = D*V' - ((final_time-t0)/2)*(Thrust.*cos(alpha)./mass - Drag./mass - g.*sin(gamma))';
ceq(3*M+1:4*M,1) = D*gamma' - ((final_time-t0)/2)*((Thrust.*sin(alpha)./(mass.*V)) - (g.*cos(gamma)./V) + (V.*cos(theta)./R))';
ceq(4*M+1:5*M,1) = D*mass' + ((final_time-t0)/2)*(Thrust./(g0*Isp))';
ceq(5*M+1,1) = Ri*Re - R(1);
ceq(5*M+2,1) = (hf/Re+1)*Re -  R(end);
ceq(5*M+3,1) = theta_i - theta(1);
ceq(5*M+4,1) = Vi*10 - V(1);
ceq(5*M+5,1) = Vf*((mu/(hf+Re))^0.5) - V(end);
ceq(5*M+6,1) = gamma_i - gamma(1);
ceq(5*M+7,1) = gamma_f - gamma(end);
ceq(5*M+8,1) = alpha_i - alpha(1);
ceq(5*M+9,1) = mass_i*m0 - mass(1);
end


