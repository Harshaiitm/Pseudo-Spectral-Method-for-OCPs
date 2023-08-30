function [c,ceq,dc,dceq] = climb_Nonlinear_func(x,N,D,mu,Thrust,r,t0,tf,Lift,Drag,g0,Isp)                        %inequality constarints
c = [];
dc = [];
dceq = [];
h = x(1:N+1);
v = x(N+2:2*N+2);
gamma = x(2*N+3:3*N+3);
mass = x(3*N+4:4*N+4);
alpha = x(4*N+5:5*N+5);
% alpha = x(4*N+6:5*N+6);


ceq = zeros(4*N+12,1);
ceq(1:N+1,1) = D*h'-((tf-t0)/2)* (v.*sin(gamma))';
ceq(N+2:2*N+2,1)=D*v' - ((tf-t0)/2)*((Thrust.*cos(alpha)-Drag)./mass - mu.*sin(gamma)./r.^2)';
ceq(2*N+3:3*N+3,1) = D*gamma' - ((tf-t0)/2)*((Thrust.*sin(alpha)+Lift)./(mass.*v)+cos(gamma).*((v./r)-mu./(v.*r.^2)))';
ceq(3*N+4:4*N+4,1) = D*mass'-((tf-t0)/2)*(Thrust./(g0.*Isp))';
ceq(4*N+5) = 0+h(1);
ceq(4*N+6) = 19995-h(end);
ceq(4*N+7) =129.314+v(1);
ceq(4*N+8) =295.092- v(end);
ceq(4*N+9) =gamma(1)-0;
ceq(4*N+10) = gamma(end)-0;
ceq(4*N+11) =19050.864-mass(1);
ceq(4*N+12) = 16000-mass(end);


end