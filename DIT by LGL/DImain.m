clc; clear all; close all;
%==============================================================================================%
N = 30; % Order of the polynomial
[nodes,weights] = LGL_nodes(N);  % Legendre_gauss_Lobatto PS Method 
D = collocD(nodes);              % Differential matrix for LGL                               
%==============================================================================================%
x = zeros(1,3*N+3);              % state and control vector assigned       
x1 = x(1:N+1);                   % position                        
x2 = x(N+2:2*N+2);               % velocity  
x3 = x(2*N+3:3*N+3);             % accleration


% initialization of state and control variables
t0 = 0;                          % initial time   
tf = 10;                         % final time
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
x1 = 5*sin(t);
x2 = 5*cos(t);
% i = linspace(1,10,N+1)
% x1 = 5*sin(i);
% x2 = 5*cos(i);
x0(1:N+1) = double(x1);
x0(N+2:2*N+2) = double(x2);
x0(2*N+3:3*N+3) = 0.1;
 

% call for Objective and Nonlinearconstraints function
Objective = Objectivee_func(x,N,weights,t);
[c, ceq, dc, dceq] = Nonlinearcon(x,N,D,t0,tf);

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:N+1) = -6;
lb(N+2:2*N+2) = -10;
lb(2*N+3:3*N+3) = -10;
ub(1:N+1) = 6;
ub(N+2:2*N+2) = 10;
ub(2*N+3:3*N+3) = 10;


%==============================================================================================%
% solver fmincon
% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000);
[x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon(x,N,D,t0,tf),options);

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);








%========================================================================================================
% Plotting 

x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
[nodes,~] = LGL_nodes(N);
t0 = 0;
tf = 10;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
% t = linspace(0,10,N+1);



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

%% polynomial
z = 8;
disp(['at time t =',num2str(z),'s']);
pointx=t';
pointy=x1;
position=lagrange(z,pointx,pointy);
disp(['Position =',num2str(position),'m']);

pointy=x2;
velocity=lagrange(z,pointx,pointy);
disp(['Velocity =',num2str(velocity),'m/s']);

pointy=x3;
acceleration=lagrange(z,pointx,pointy);
disp(['Acceleration =',num2str(acceleration),'m/s^2']);







