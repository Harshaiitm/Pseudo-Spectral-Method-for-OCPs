function [c,ceq,dc,dceq] = climb_Nonlinear_func_LGR(x,M,D,problem)                        %inequality constarints
c = [];
dc = [];
dceq = [];

h = x(1:M);
v = x(M+1:2*M);
gamma = x(2*M+1:3*M);
mass = x(3*M+1:4*M);
alpha = x(4*M+1:5*M);
final_time = x(5*M+1);

Re = problem.Re;
mu = problem.mu;
S = problem.S;
g0 = problem.g0;
Isp = problem.Isp;
t0 = problem.t0;
x0 = problem.x0;
hi = problem.hi;
hf = problem.hf;
vi = problem.vi;
vf = problem.vf;
gamma_i = problem.gamma_i;
gamma_f = problem.gamma_f;
mass_i = problem.mass_i;

r = h + Re;
[rho,sos]=atm_data(h);
Mach =  v./sos;


[Clalpha,CD0,eta] = aero_data(Mach);
Thrust = thrust_available(h,Mach);


CL = Clalpha.*alpha;
CD = CD0 + eta.*Clalpha.*alpha.^2;
q = 0.5.*rho.*v.*v;
Drag = q.*S.*CD;
Lift = q.*S.*CL;



ceq = zeros(4*M+3,1);
ceq(1:M,1) = D*[hi h]'-((final_time-t0)/2)* (v.*sin(gamma))';
ceq(M+1:2*M,1)= D*[vi v]' - ((final_time-t0)/2)*((Thrust.*cos(alpha)-Drag)./mass - mu.*sin(gamma)./r.^2)';
ceq(2*M+1:3*M,1) = D*[gamma_i gamma]' - ((final_time-t0)/2)*((Thrust.*sin(alpha)+Lift)./(mass.*v)+cos(gamma).*((v./r)-mu./(v.*r.^2)))';
ceq(3*M+1:4*M,1) = D*[mass_i mass]'+((final_time-t0)/2)*(Thrust./(g0.*Isp))';
ceq(4*M+1) = hf - h(M);
ceq(4*M+2) = vf - v(M);
ceq(4*M+3) = gamma_f - gamma(M);

end