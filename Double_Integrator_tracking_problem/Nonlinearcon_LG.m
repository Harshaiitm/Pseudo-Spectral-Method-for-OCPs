function [c, ceq, dc, dceq] = Nonlinearcon_LG(x,N,D,t0,tf,weights)
c = [];
dc = [];
dceq = [];
x1 = x(1:N);
x1(N+1)=x1(1)+weights(1:N)'*(x1(1:N))';
x2 = x(N+2:2*N+1);
x2(N+1)=x2(1)+weights(1:N)'*(x2(1:N))';
x3 = x(2*N+3:3*N+2);
x3(N+1)=x3(1)+weights(1:N)'*(x3(1:N))';

ceq = zeros(2*N+5,1);
ceq(1:N,1) = D * x1' - ((tf-t0)/2)*(x2(1:N))';
ceq(N+1,1) =x1(N+1)-x1(1)-weights(1:N)'*(x1(1:N))';
ceq(N+2,1) =x2(N+1)-x2(1)-weights(1:N)'*(x2(1:N))';
ceq(N+3:2*N+2,1) = D * x2' -((tf-t0)/2)*(x3(1:N))';
ceq(2*N+3,1) =x3(N+1)-x3(1)-weights(1:N)'*(x3(1:N))';
ceq(2*N+4,1) = x(1)-0;            % starting point
ceq(2*N+5,1) = x(N+2)-5;
end
