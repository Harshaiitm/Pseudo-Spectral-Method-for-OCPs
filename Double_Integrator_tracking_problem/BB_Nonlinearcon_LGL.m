function [c, ceq, dc, dceq] = BB_Nonlinearcon_LGL(x,M,D,problem)

xi = problem.xi;
xf = problem.xf;
vi = problem.vi;
vf = problem.vf;
x0 = problem.x0;
t0 = problem.t0;

c = [];
dc = [];
dceq = [];
x1 = x(1:M);                        % position vector
x2 = x(M+1:2*M);                    % velocity vector
x3 = x(2*M+1:3*M);                  % Acceleration vector
x4 = x(3*M+1);                      % final time


ceq = zeros(2*M+4,1);
ceq(1:M,1) = D * x1' - ((x4-t0)/2)*x2';
ceq(M+1:2*M,1) = D * x2' -((x4-t0)/2)*x3';
ceq(2*M+1) = x1(1) - xi;
ceq(2*M+2) = x1(M) - xf;
ceq(2*M+3) = x2(1) - vi;
ceq(2*M+4) = x2(M) - vf;

end
