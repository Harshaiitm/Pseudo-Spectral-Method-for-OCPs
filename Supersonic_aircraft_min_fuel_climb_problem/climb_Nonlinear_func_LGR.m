function [c,ceq,dc,dceq] = climb_Nonlinear_func_LGR(x,N,D,problem)                        %inequality constarints
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



ceq = zeros(4*N+7,1);
ceq(1:N,1) = D*h'-((final_time-t0)/2)* (v(1:N).*sin(gamma(1:N)))';
ceq(N+1:2*N,1)=D*v' - ((final_time-t0)/2)*((Thrust(1:N).*cos(alpha(1:N))-Drag(1:N))./mass(1:N) - mu.*sin(gamma(1:N))./r(1:N).^2)';
ceq(2*N+1:3*N,1) = D*gamma' - ((final_time-t0)/2)*((Thrust(1:N).*sin(alpha(1:N))+Lift(1:N))./(mass(1:N).*v(1:N))+cos(gamma(1:N)).*((v(1:N)./r(1:N))-mu./(v(1:N).*r(1:N).^2)))';
ceq(3*N+1:4*N,1) = D*mass'+((final_time-t0)/2)*(Thrust(1:N)./(g0.*Isp))';
ceq(4*N+1) = 0-h(1);
ceq(4*N+2) = h(N)-19995;
ceq(4*N+3) =129.314-v(1);
ceq(4*N+4) =v(N)-295.092;
ceq(4*N+5) =0-gamma(1);
ceq(4*N+6) = gamma(N)-0;
ceq(4*N+7) =19050.864-mass(1);



end