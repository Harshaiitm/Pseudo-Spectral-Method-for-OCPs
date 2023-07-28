clc; clear all; close all;
N = 20; % Order of the polynomial
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
disp(x(1:N))

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




