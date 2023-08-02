clc; clear all; close all;
N = 20; % Order of the polynomial

x = zeros(1,3*N+4); 
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);

[nodes,weights,~] = LGR_nodes(N);
size(nodes);
nodes = flip(-nodes);
nodes(1) = -1;
weights = flip(weights);
D = collocD(nodes);
D(N+1,:) = [];

t0 = 0;
tf = 35;


Objective = BB_Objective_func_LGR(x,N);
[c, ceq, dc, dceq] = BB_Nonlinearcon_func_LGR(x,N,D,t0,tf);



x0(1) = 0;
i =linspace(0,300,N-1);
x0(2:N) = double(i);
x0(N+1) = 300;
x0(N+2) = 0;
j =linspace(0,200,N-1);
x0(N+3:2*N+1) = double(j);
x0(2*N+2) = 0;
x0(2*N+3:3*N+3) = 0.5;
x0(3*N+4) = 20;

A = [];
b = [];
Aeq = [];
beq = [];
lb = x;
lb(1:N+1) = -10;
lb(N+2:2*N+2) = -200;
lb(2*N+3:3*N+3) = -2;
lb(3*N+4) = 0;
ub = x;
ub(1:N+1) = 300;
ub(N+2:2*N+2) = 200;
ub(2*N+3:3*N+3) = 1;
ub(3*N+4) = 35;


% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000);
[x,fval,ef,output] = fmincon(@BB_Objective_func_LGR,x0,A,b,Aeq,beq,lb,ub,@(x)BB_Nonlinearcon_func_LGR(x,N,D,t0,tf),options);

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);



%========================================================================================================
% 
% Plotting 

x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
% t = linspace(0,40,N+1);



figure(1)
plot(t, x1);
hold on
plot(t, x2);
xlabel('Time (s)');
ylabel('State Variables');
title('Double Integrator Tracking Problem');
legend({'Positon(x1)','Velocity(x2)'},Location="northeast");
grid on
hold off


figure(2)
plot(t,x3);
xlabel('Time (s)');
ylabel('countrol variables (N)');
legend({'control variable'},Location="northeast");
title('Double integrator tracking problem');

%% Polynomial
z_value = 35;
disp(['at time t =',num2str(z_value),'s']);
coeff_P = polyfit(t,x1,N);
sympref('FloatingPointOutput',true);
syms z
polynomial_p = sym(0);
for i= 1:N+1
    polynomial_p =vpa(polynomial_p+coeff_P(i)*z^(N-i+1));
end
% polynomial = polyval(coeff', z);
disp(['Position Equation:', char(polynomial_p)]);
position = subs(polynomial_p,z,z_value);
disp(['Position:', char(position),'m']);

coeff_V = polyfit(t,x2,N);
polynomial_V = sym(0);
for i= 1:N+1
    polynomial_V =vpa(polynomial_V+coeff_V(i)*z^(N-i+1));
end
disp(['Velocity Equation:', char(polynomial_V)]);
velocity = subs(polynomial_V,z,z_value);
disp(['velocity:', char(velocity),'m/s']);

coeff_A = polyfit(t,x3,N);
polynomial_A = sym(0);
for i= 1:N+1
    polynomial_A =vpa(polynomial_A+coeff_A(i)*z^(N-i+1));
end
disp(['Acceleration Equation:', char(polynomial_A)]);
accleration = subs(polynomial_A,z,z_value);
disp(['Acceleration:', char(accleration),'m/s^2']);



