clc; clear all; close all;
N = 20; % Order of the polynomial
[nodes,~] = LGL_nodes(N);
D = collocD(nodes);

x = zeros(1,3*N+4); 
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);
t0 = 0;
tf = 35;


Objective = Objectivee_func(x,N);
[c, ceq, dc, dceq] = Nonlinearcon_func(x,N,D,t0,tf);



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
[x,fval,ef,output] = fmincon(@Objectivee_func,x0,A,b,Aeq,beq,lb,ub,@(x)Nonlinearcon_func(x,N,D,t0,tf),options);

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

% Lagrange interpolation
z_value = 8;  % at time in seconds
disp(['at time t =',num2str(z_value),'s']);
pointx=t';
pointy=x1;
syms z
Position_equation = lagrange_interpolation(pointx,pointy);
disp(['Position Equation:', char(Position_equation)]);
position = subs(Position_equation,z,z_value);
disp(['Position =',char(position),'m']);

pointy=x2;
velocity_equation=lagrange_interpolation(pointx,pointy);
disp(['Velocity Equation:', char(velocity_equation)]);
velocity = subs(velocity_equation,z,z_value);
disp(['Velocity =',char(velocity),'m/s']);

pointy=x3;
acceleration_equation=lagrange_interpolation(pointx,pointy);
disp(['Acceleration =',char(acceleration_equation),'m/s^2']);
acceleration = subs(acceleration_equation,z,z_value);
disp(['acceleration =',char(acceleration),'m/s']);


