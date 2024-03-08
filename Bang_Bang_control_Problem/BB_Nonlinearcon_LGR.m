function [c, ceq, dc, dceq] = BB_Nonlinearcon_LGR(x,x0,M,D,t0)

dc = [];
dceq = [];
x1 = x(1:M);                        % position vector
x2 = x(M+1:2*M);                    % velocity vector
x3 = x(2*M+1:3*M);                  % Acceleration vector
x4 = x(3*M+1);                      % final time

c = zeros(M,1);
c(1:M,1) = 0-x1;

ceq = zeros(2*M+3,1);
ceq(1:M,1) = D*[x0(1) x1]' - ((x4-t0)/2) * x2';
ceq(M+1:2*M,1) = D*[x0(M+1) x2]' -((x4-t0)/2) * x3';
ceq(2*M+1) = x1(end)-300;
ceq(2*M+2) = x2(end)-0;

end
