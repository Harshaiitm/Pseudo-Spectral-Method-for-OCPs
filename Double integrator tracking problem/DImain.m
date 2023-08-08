%--------------------------------------------------------------------------
% DImain.m
% - Main script for the Double integrator Tracking problem
% - Define the mesh and method, then solve and plot solution
%--------------------------------------------------------------------------
% DImain
%--------------------------------------------------------------------------
clc; clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGR';   % either LGL or LG or LGR
N = 40;     % Order of the polynomial

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
    if strcmp(PS_method,'LGL')
       Objective = Objectivee_func_LGL(x,N,weights,t);
       [c, ceq, dc, dceq] = Nonlinearcon_LGL(x,N,D,t0,tf);
       [x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGL(x,N,D,t0,tf),options);
    elseif strcmp(PS_method,'LG')
       Objective = Objectivee_func_LG(x,N,weights,t); 
       [c, ceq, dc, dceq] = Nonlinearcon_LG(x,N,D,t0,tf);
       [x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LG(x,N,D,t0,tf),options);
    elseif strcmp(PS_method,'LGR')
       Objective = Objectivee_func_LGR(x,N,weights,t);
       [c, ceq, dc, dceq] = Nonlinearcon_LGR(x,N,D,t0,tf);
       [x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGR(x,N,D,t0,tf),options);
    end 

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);



%========================================================================================================
% Plotting 

x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
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







