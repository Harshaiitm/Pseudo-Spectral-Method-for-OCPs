function [c, ceq, dc, dceq] = BB_Nonlinearcon_LGR(x,M,D,problem)

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


ceq = zeros(2*M+2,1);
ceq(1:M,1) = D*[xi x1]' - ((x4-t0)/2) * x2';
ceq(M+1:2*M,1) = D*[vi x2]' -((x4-t0)/2) * x3';
ceq(2*M+1) = xf - x1(M);
ceq(2*M+2) = vf - x2(M);

end
