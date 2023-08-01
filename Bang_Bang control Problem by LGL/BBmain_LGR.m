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




