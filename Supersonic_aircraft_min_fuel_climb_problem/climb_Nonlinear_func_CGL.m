function [c,ceq,dc,dceq] = climb_Nonlinear_func_CGL(x,N,D,problem)                        %inequality constarints
c = [];
dc = [];
dceq = [];

h = x(1:N+1);
v = x(N+2:2*N+2);
gamma = x(2*N+3:3*N+3);
final_time = x(3*N+4);
mass = x(3*N+5:4*N+5);
alpha = x(4*N+6:5*N+6);

Re = problem.Re;
mu = problem.mu;
S = problem.S;
g0 = problem.g0;
Isp = problem.Isp;
t0 = problem.t0;

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



ceq = zeros(4*N+11,1);
ceq(1:N+1,1) = D*h'-((final_time-t0)/2)* (v.*sin(gamma))';
ceq(N+2:2*N+2,1)=D*v' - ((final_time-t0)/2)*((Thrust.*cos(alpha)-Drag)./mass - mu.*sin(gamma)./r.^2)';
ceq(2*N+3:3*N+3,1) = D*gamma' - ((final_time-t0)/2)*((Thrust.*sin(alpha)+Lift)./(mass.*v)+cos(gamma).*((v./r)-mu./(v.*r.^2)))';
ceq(3*N+4:4*N+4,1) = D*mass'+((final_time-t0)/2)*(Thrust./(g0.*Isp))';
ceq(4*N+5) = 0-h(1);
ceq(4*N+6) = h(end)-19995;
ceq(4*N+7) =129.314-v(1);
ceq(4*N+8) =v(end)-295.092;
ceq(4*N+9) =0-gamma(1);
ceq(4*N+10) = gamma(end)-0;
ceq(4*N+11) =19050.864-mass(1);


end