function [c, ceq, dc, dceq] = Nonlinearcon_LG(x,N,D,t0,tf,weights)

c = [];
dc = [];
dceq = [];
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);

c =zeros(N+1,1);
c(1:N,1) = 0-x1(1:N);
c(N+1,1) = 0-x1(N+1);


ceq = zeros(2*N+4,1);
ceq(1:N,1) = D * x1' - ((tf-t0)/2)*(x2(1:N))';
% ceq(N+1,1) =x1(N+1)-x2(1)-weights(1:N)'*(x1(1:N))';
ceq(N+1:2*N,1) = D * x2' -((tf-t0)/2)*(x3(1:N))';
% ceq(2*N+2,1) =x3(N+1)-x3(1)-weights(1:N)'*(x2(1:N))';
ceq(2*N+1) = x1(end)-300;
ceq(2*N+2) = x2(1)-0;
ceq(2*N+3) = x2(end)-0;
ceq(2*N+4) =x1(end)-(x3(end)/2)*(x4^2);
end
