% Solved using an OPTI NLP solver.
clc
fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;    
lb = [-1.5;-3];
ub = [4;3];
x0 = [0;0];


opts = optiset('solver','IPOPT');
[x,fval,exitflag,info] = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts);



