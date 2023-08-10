%--------------------------------------------------------------------------
% BBmain.m
% - Main script for the Bang_Bang control problem
% - Define the PS method and Order of the polynomial, then solve and plot solution
%--------------------------------------------------------------------------
% BBmain
%--------------------------------------------------------------------------
clc; clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N = 50;     % Order of the polynomial
addpath('C:\Users\Harshad\OneDrive\Desktop\DIT\PS_methods')  % add the PS_method file directory

     if  strcmp(PS_method,'LGL')
        [nodes,weights] = LGL_nodes(N); % calculate scaled node locations and weights
        D=collocD(nodes); % segment differentiation matrix
    elseif strcmp(PS_method,'LG')
        nodes(1) =-1;
        [nodes(2:N+1,1),weights]=LG_nodes(N,-1,1); % calculate scaled node locations and weights
        nodes(2:N+1,1) = flip(nodes(2:N+1,1)); 
        D=collocD(nodes); % segment differentiation matrix
        D(N+1,:) = [];
    elseif strcmp(PS_method,'LGR')
        [nodes,weights] = LGR_nodes(N); % calculate scaled node locations and weights
        D = collocD(nodes); % segment differentiation matrix
        nodes = flip(-nodes);
        nodes(1) = -1;
        weights = flip(weights);
        D(N+1,:) = [];
    end 

%==============================================================================================%

x = zeros(1,3*N+4); 
x1 = x(1:N+1);              % position vector
x2 = x(N+2:2*N+2);          % velocity vector
x3 = x(2*N+3:3*N+3);        % Acceleration vector
x4 = x(3*N+4);              % final time
t0 = 0;
tf = 35;
  


% initialization of state and control variables
x0(1) = 0;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
x0(1:N+1) = t*8.5714;
% x0(2:N) =linspace(0,300,N-1);
% x0(N+1) = 300;
x0(N+2) = 0;
j =linspace(0,200,N-1);
% x0(N+3:2*N+1) = double(j);
x0(2:N) = t(2:N)*1.5;
x0(2*N+2) = 0;
x0(2*N+3:3*N+3) =2.5;
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
    if strcmp(PS_method,'LGL')
       [c, ceq, dc, dceq] = Nonlinearcon_LGL(x,N,D,t0,tf);
       [x,fval,ef,output] = fmincon(@(x) Objective_func(x),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGL(x,N,D,t0,tf),options);
    elseif strcmp(PS_method,'LG') 
       [c, ceq, dc, dceq] = Nonlinearcon_LG(x,N,D,t0,tf);
       [x,fval,ef,output] = fmincon(@(x) Objective_func(x),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LG(x,N,D,t0,tf),options);
    elseif strcmp(PS_method,'LGR')
       [c, ceq, dc, dceq] = Nonlinearcon_LGR(x,N,D,t0,tf);
       [x,fval,ef,output] = fmincon(@(x) Objective_func(x),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGR(x,N,D,t0,tf),options);
    end 

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

% Lagrange interpolation
z_value = t;  % at time in seconds
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





figure(1)
plot(t, position);
hold on
plot(t, velocity);
xlabel('Time (s)');
ylabel('State Variables');
title('Double Integrator Tracking Problem');
legend({'Positon(x1)','Velocity(x2)'},Location="northeast");
grid on
hold off


figure(2)
plot(t,acceleration);
xlabel('Time (s)');
ylabel('countrol variables (N)');
legend({'control variable'},Location="northeast");
title('Double integrator tracking problem');




