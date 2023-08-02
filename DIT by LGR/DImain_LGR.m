clc; clear all; close all;
N = 40; % Order of the polynomial
x = zeros(1,3*N+3); 
x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
[nodes,weights,~] = LGR_nodes(N);
size(nodes);
nodes = flip(-nodes);
nodes(1) = -1;
weights = flip(weights);
D = collocD(nodes);
D(N+1,:) = [];
t0 =0;
tf =10;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;

Objective = Objectivee_func_LGR(x,N,weights,t);
[c, ceq, dc, dceq] = Nonlinearcon_LGR(x,N,D,t0,tf);


x1 = 5*sin(t);
x2 = 5*cos(t);
% i = linspace(1,10,N+1)
% x1 = 5*sin(i);
% x2 = 5*cos(i);
x0(1:N+1) = double(x1);
x0(N+2:2*N+2) = double(x2);
x0(2*N+3:3*N+3) = 0.1;

A = [];
b = [];
Aeq = [];
beq = [];

lb(1:N+1) = -6;
lb(N+2:2*N+2) = -10;
lb(2*N+3:3*N+3) = -10;
ub(1:N+1) = 6;
ub(N+2:2*N+2) = 10;
ub(2*N+3:3*N+3) = 10;

% Start the timer
tic;

options =  optimoptions ('fmincon','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000, 'Algorithm','sqp-legacy' );
[x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x)Nonlinearcon_LGR(x,N,D,t0,tf),options);

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);


%========================================================================================================
% 
% Plotting 

x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);

% t = linspace(0,10,N+1);



figure(1)
plot(t, x1);
hold on
plot(t, x2);
xlabel('Time (s)');
ylabel('State Variables');
title('Double Integrator Tracking Problem by LGR PS Method');
legend({'Positon(x1)','Velocity(x2)'},Location="northeast");
hold off
  

figure(2)
plot(t,x3);
xlabel('Time (s)');
ylabel('countrol variables (N)');
legend({'control variable'},Location="northeast");
title('Double integrator tracking problem by LGR PS Method');

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

