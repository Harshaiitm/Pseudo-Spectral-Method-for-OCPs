function [c,ceq,dc,dceq] = Normal_Tangential_Polar_Rocket_Nonlinear_func_LGL(x,N,D,problem)   %inequality constarints                      
% c = [];
dc = [];
dceq = [];

% Decision veriables
R = x(1:N+1);               
theta = x(N+2:2*N+2);
V = x(2*N+3:3*N+3);
gamma = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust = x(5*N+6:6*N+6);
alpha = x(6*N+7:7*N+7);
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
g = mu./(R).^2;
q = 0.5*rho.*(V.^2);
Drag = q.* A_ref *CD;
a_sen_v = (Thrust.* cos(alpha) - Drag)./mass;
a_sen_gamma = (Thrust.* sin(alpha))./mass;
a_sen_mag = (a_sen_v.^2 + a_sen_gamma.^2).^0.5;

c = zeros(3*N+3,1);
c(1:N+1,1) = Thrust.^2 - Thrust_max^2;
c(N+2:2*N+2,1) = a_sen_mag.^2 - a_sen_max^2;
c(2*N+3:3*N+3,1) = q - q_max;

ceq = zeros(5*N+14,1);
ceq(1:N+1,1) = D*R' - ((final_time-t0)/2)*(V.*sin(gamma))';
ceq(N+2:2*N+2,1) = D*theta' - ((final_time-t0)/2)*(V.*cos(gamma)./R)';
ceq(2*N+3:3*N+3,1) = D*V' - ((final_time-t0)/2)*(Thrust.*cos(alpha)./mass - Drag./mass - g.*sin(gamma))';
ceq(3*N+4:4*N+4,1) = D*gamma' - ((final_time-t0)/2)*((Thrust.*sin(alpha)./(mass.*V) + g.*cos(gamma)./V) + V.*cos(theta)./R)';
ceq(4*N+5:5*N+5,1) = D*mass' + ((final_time-t0)/2)*(Thrust./(g0*Isp))';
ceq(5*N+6,1) = Re - R(1);
ceq(5*N+7,1) = (hf+Re) -  R(end);
ceq(5*N+8,1) = 0 - theta(1);
ceq(5*N+9,1) = 0 - V(1);
ceq(5*N+10,1) = (mu/(hf+Re))^0.5- V(end);
ceq(5*N+11,1) = 0 - gamma(1);
ceq(5*N+12,1) = 0 - gamma(end);
ceq(5*N+13,1) = 0 - alpha(1);
ceq(5*N+14,1) = m0 - mass(1);
end


