clc; clear all; close all;
N = 100; % Order of the polynomial
nodes(1) =-1;
% nodes(N+2) =1;
[nodes(2:N+1,1),weights]=LG_nodes(N,-1,1);
nodes(2:N+1,1) = flip(nodes(2:N+1,1));
D = collocD(nodes);
D(N+1,:) = []; 
t0 =0;
tf =10;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
x = zeros(1,3*N+3); 
x1 = x(1:N+1);
% x1(N+2) = x1(1)+((tf-t0)/2)*(x1(1:N))*weights;
x2 = x(N+2:2*N+2);
% x2(N+2) = x2(1)+((tf-t0)/2)*(x2(1:N))*weights;
x3 = x(2*N+3:3*N+3);
% x3(N+2) = x3(1)+((tf-t0)/2)*(x3(1:N))*weights;

Objective = Objectivee_func_LG(x,N,weights,t,t0,tf);
[c, ceq, dc, dceq] = Nonlinearcon_LG(x,N,D,t0,tf,weights);


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
[x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x)Nonlinearcon_LG(x,N,D,t0,tf,weights),options);


% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);


%========================================================================================================
% 
% Plotting 

x1 = x(1:N+1);
x1(N+2) = x1(1)+((tf-t0)/2)*(x1(1:N))*weights;
x2 = x(N+2:2*N+2);
x2(N+2) = x2(1)+((tf-t0)/2)*(x2(1:N))*weights;
x3 = x(2*N+3:3*N+3);
x3(N+2) = x3(1)+((tf-t0)/2)*(x3(1:N))*weights;

nodes(N+2) = 1;
t(N+2) = ((tf-t0)/2).*nodes(N+2)+(tf+t0)/2;

% t = linspace(0,10,N+1);



figure(1)
plot(t, x1);
hold on
plot(t, x2);
xlabel('Time (s)');
ylabel('State Variables');
title('Double Integrator Tracking Problem by LG PS Method');
legend({'Positon(x1)','Velocity(x2)'},Location="northeast");
hold off
  
 
figure(2)
plot(t,x3);
xlabel('Time (s)');
ylabel('countrol variables (N)');
legend({'control variable'},Location="northeast");
title('Double integrator tracking problem by LG PS Method');




